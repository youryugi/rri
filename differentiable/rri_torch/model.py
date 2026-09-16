"""Assembles slope routing, river routing, river-slope exchange,
infiltration and evaporation into one time-stepping model, preserving the
operator-splitting order used by the main loop in RRI.f90:

  1. river routing over the full step dt (RK, ignoring slope exchange)
  2. slope routing over the full step dt (RK, rainfall forcing only)
  3. evapotranspiration (explicit, post-routing)
  4. river <-> slope exchange (weir formulas)
  5. infiltration (Green-Ampt, explicit, post-exchange)
  6. domain-boundary ("sink") cells drain to zero

Not implemented (see README for the full list and rationale): adaptive
RKF45 time-stepping, groundwater, dams, diversions, levee breaks, sediment
transport, boundary-condition files, and non-rectangular cross-sections.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import torch

from .geometry import Grid
from .slope import SlopeParams, slope_rhs
from .river import RiverParams, river_rhs
from .exchange import ExchangeGeometry, river_slope_exchange
from .processes import infiltration, evaporation
from .integrate import integrate_fixed


@dataclass
class RRIParams:
    slope: SlopeParams
    river: RiverParams
    ksv: torch.Tensor             # (ny, nx) Green-Ampt saturated hydraulic conductivity [m/s]
    faif: torch.Tensor            # (ny, nx) Green-Ampt wetting-front suction [m]
    infilt_limit: torch.Tensor    # (ny, nx) cap on cumulative infiltration [m], < 0 = unlimited
    evp_switch: int = 0           # 0 = no ET, 1 = ET can draw from gampt_ff, 2 = ET from hs only


@dataclass
class StepDiagnostics:
    qr: torch.Tensor          # (n_riv,) discharge to downstream cell at end of step [m^3/s]
    qrs: torch.Tensor         # (ny, nx) slope->river exchange rate [m/s]
    aevp: torch.Tensor        # (ny, nx) actual evapotranspiration rate [m/s]
    infilt_rate: torch.Tensor  # (ny, nx) infiltration rate [m/s]
    drained_volume: torch.Tensor  # scalar, water removed at sink cells this step [m^3]


class RRIModel:
    def __init__(self, grid: Grid, n_substeps_slope: int = 4, n_substeps_river: int = 4):
        self.grid = grid
        self.n_substeps_slope = n_substeps_slope
        self.n_substeps_river = n_substeps_river

    def initial_state(self, dtype: torch.dtype = torch.float64):
        g = self.grid
        hs = torch.zeros((g.ny, g.nx), dtype=dtype)
        gampt_ff = torch.zeros((g.ny, g.nx), dtype=dtype)
        hr = torch.zeros((g.n_riv,), dtype=dtype)
        return hs, gampt_ff, hr

    def step(
        self,
        hs: torch.Tensor,
        gampt_ff: torch.Tensor,
        hr: torch.Tensor,
        params: RRIParams,
        rain: torch.Tensor,
        pet: Optional[torch.Tensor],
        dt: float,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, StepDiagnostics]:
        g = self.grid

        # 1. river routing (independent of slope for the duration of dt)
        def river_rhs_fn(hr_):
            dhrdt, _ = river_rhs(hr_, g, params.river)
            return dhrdt

        hr1 = integrate_fixed(river_rhs_fn, hr, dt, self.n_substeps_river)
        _, qr_end = river_rhs(hr1, g, params.river)

        # 2. slope routing (rainfall forcing only)
        def slope_rhs_fn(hs_):
            dhsdt, _ = slope_rhs(hs_, g.zb, g.domain, params.slope, rain, g.dx, g.dy, g.area)
            return dhsdt

        hs1 = integrate_fixed(slope_rhs_fn, hs, dt, self.n_substeps_slope)

        # 3. evapotranspiration
        if params.evp_switch and pet is not None:
            hs1, gampt_ff1, aevp = evaporation(hs1, gampt_ff, pet, dt, params.evp_switch)
        else:
            gampt_ff1 = gampt_ff
            aevp = torch.zeros_like(hs1)

        # 4. river <-> slope exchange
        exch_geom = ExchangeGeometry(
            riv_mask=g.riv_mask,
            depth=g.scatter_from_riv(params.river.depth),
            height=g.scatter_from_riv(params.river.height),
            width=g.scatter_from_riv(params.river.width),
            len_riv=g.scatter_from_riv(g.len_riv),
            area=g.area,
        )
        hr2_grid_in = g.scatter_from_riv(hr1)
        hr2_grid, hs2, qrs = river_slope_exchange(hr2_grid_in, hs1, exch_geom, dt)
        hr2 = g.gather_to_riv(hr2_grid)

        # 5. infiltration
        hs3, gampt_ff2, infilt_rate = infiltration(
            hs2, gampt_ff1, params.ksv, params.faif, params.slope.gammaa, params.infilt_limit, dt
        )

        # 6. drain sink (open-boundary / outlet) cells
        sink_riv = g.gather_to_riv(g.sink)
        channel_area = torch.clamp(params.river.width * g.len_riv, min=1e-12)
        drained = (hs3[g.sink] * g.area).sum() + (hr2[sink_riv] * channel_area[sink_riv]).sum()
        hs4 = torch.where(g.sink, torch.zeros_like(hs3), hs3)
        hr3 = torch.where(sink_riv, torch.zeros_like(hr2), hr2)

        diag = StepDiagnostics(qr=qr_end, qrs=qrs, aevp=aevp, infilt_rate=infilt_rate, drained_volume=drained)
        return hs4, gampt_ff2, hr3, diag

    def simulate(
        self,
        rain_seq: torch.Tensor,
        params: RRIParams,
        dt: float,
        pet_seq: Optional[torch.Tensor] = None,
        outlet_riv_index: Optional[int] = None,
        init_state: Optional[tuple] = None,
        record_full_state: bool = False,
    ) -> dict:
        """Run the model over a forcing time series.

        rain_seq: (T, ny, nx) rainfall rate [m/s] for each step.
        pet_seq: optional (T, ny, nx) potential ET rate [m/s].
        outlet_riv_index: river-cell index to report a discharge series
            for (typically the catchment outlet); if None, no `qr_outlet`
            key is returned.
        record_full_state: if True, also stack the full `hs`/`hr` history
            (memory-heavy for long runs / big grids).

        Returns a dict of stacked tensors plus the final state under
        "hs_final", "gampt_ff_final", "hr_final".
        """
        hs, gampt_ff, hr = init_state if init_state is not None else self.initial_state(dtype=rain_seq.dtype)

        qr_outlet_series = []
        drained_series = []
        hs_series = []
        hr_series = []

        T = rain_seq.shape[0]
        for t in range(T):
            pet_t = pet_seq[t] if pet_seq is not None else None
            hs, gampt_ff, hr, diag = self.step(hs, gampt_ff, hr, params, rain_seq[t], pet_t, dt)
            drained_series.append(diag.drained_volume)
            if outlet_riv_index is not None:
                qr_outlet_series.append(diag.qr[outlet_riv_index])
            if record_full_state:
                hs_series.append(hs)
                hr_series.append(hr)

        out = {
            "hs_final": hs,
            "gampt_ff_final": gampt_ff,
            "hr_final": hr,
            "drained_volume": torch.stack(drained_series),
        }
        if outlet_riv_index is not None:
            out["qr_outlet"] = torch.stack(qr_outlet_series)
        if record_full_state:
            out["hs"] = torch.stack(hs_series)
            out["hr"] = torch.stack(hr_series)
        return out
