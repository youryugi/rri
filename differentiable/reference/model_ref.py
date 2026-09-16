"""Assembles the reference port into the same operator-split time step as
RRI.f90's main loop (see rri_torch/model.py for the same list applied to
the differentiable version):

  1. river routing over dt      (RK, ignoring slope exchange)
  2. slope routing over dt      (RK, rainfall forcing only)
  3. evapotranspiration
  4. river <-> slope exchange
  5. infiltration
  6. drain sink cells

Supports two time-stepping modes:
  - "fixed": fixed-step RK4 with `n_substeps_*`, directly comparable to
    `rri_torch.model.RRIModel` at the same (dt, n_substeps).
  - "adaptive": the real embedded RKF45 adaptive stepping RRI.f90 uses,
    the closest thing to "what the actual Fortran solver would compute"
    this reference can offer without the compiled binary + a CGNS input
    (see the differentiable/README's "Validating against the Fortran
    model" section for why that path isn't available here).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .geometry_ref import RefGrid
from .physics_ref import (
    slope_rhs_ref, river_rhs_ref, river_slope_exchange_ref,
    infiltration_ref, evaporation_ref,
)
from .integrate_ref import (
    integrate_fixed_grid, integrate_adaptive_grid,
    integrate_fixed_dict, integrate_adaptive_dict,
)


@dataclass
class RefParams:
    slope: dict            # {"ns","ka","beta","da","dm","soildepth","gammaa","min_wc4latflow"} -> (ny,nx) arrays
    river: dict             # {"width","depth","height"} -> {(i,j): value}, {"ns_river","outlet_slope"} -> float
    ksv: np.ndarray
    faif: np.ndarray
    infilt_limit: np.ndarray
    evp_switch: int = 0


@dataclass
class RefDiagnostics:
    qr: dict
    qrs: dict
    aevp: np.ndarray
    infilt_rate: np.ndarray
    drained_volume: float


class RefModel:
    def __init__(
        self,
        grid: RefGrid,
        mode: str = "fixed",
        n_substeps_slope: int = 10,
        n_substeps_river: int = 10,
        eps: float = 0.01,
        ddt_min_riv: float = 0.1,
        ddt_min_slo: float = 1.0,
    ):
        assert mode in ("fixed", "adaptive")
        self.grid = grid
        self.mode = mode
        self.n_substeps_slope = n_substeps_slope
        self.n_substeps_river = n_substeps_river
        self.eps = eps
        self.ddt_min_riv = ddt_min_riv
        self.ddt_min_slo = ddt_min_slo

    def initial_state(self):
        g = self.grid
        hs = np.zeros((g.ny, g.nx))
        gampt_ff = np.zeros((g.ny, g.nx))
        hr = {k: 0.0 for k in g.down_riv}
        return hs, gampt_ff, hr

    def _integrate_river(self, rhs, hr, dt):
        if self.mode == "fixed":
            return integrate_fixed_dict(rhs, hr, dt, self.n_substeps_river)
        return integrate_adaptive_dict(rhs, hr, dt, self.eps, self.ddt_min_riv)

    def _integrate_slope(self, rhs, hs, dt):
        if self.mode == "fixed":
            return integrate_fixed_grid(rhs, hs, dt, self.n_substeps_slope)
        return integrate_adaptive_grid(rhs, hs, dt, self.grid.domain, self.eps, self.ddt_min_slo)

    def step(self, hs, gampt_ff, hr, params: RefParams, rain, pet, dt):
        g = self.grid

        def river_rhs_fn(hr_):
            dhrdt, _ = river_rhs_ref(hr_, g, params.river)
            return dhrdt

        hr1 = self._integrate_river(river_rhs_fn, hr, dt)
        _, qr_end = river_rhs_ref(hr1, g, params.river)

        def slope_rhs_fn(hs_):
            return slope_rhs_ref(hs_, g, params.slope, rain)

        hs1 = self._integrate_slope(slope_rhs_fn, hs, dt)

        if params.evp_switch and pet is not None:
            hs2 = np.zeros_like(hs1)
            gampt_ff1 = np.zeros_like(gampt_ff)
            aevp = np.zeros_like(hs1)
            for i in range(g.ny):
                for j in range(g.nx):
                    if not g.domain[i, j]:
                        continue
                    hs2[i, j], gampt_ff1[i, j], aevp[i, j] = evaporation_ref(
                        hs1[i, j], gampt_ff[i, j], pet[i, j], dt, params.evp_switch
                    )
        else:
            hs2 = hs1
            gampt_ff1 = gampt_ff
            aevp = np.zeros_like(hs1)

        hr2, hs3, qrs = river_slope_exchange_ref(hr1, hs2, g, params.river, dt)

        hs4 = np.zeros_like(hs3)
        gampt_ff2 = np.zeros_like(gampt_ff1)
        infilt_rate = np.zeros_like(hs3)
        for i in range(g.ny):
            for j in range(g.nx):
                if not g.domain[i, j]:
                    continue
                hs4[i, j], gampt_ff2[i, j], infilt_rate[i, j] = infiltration_ref(
                    hs3[i, j], gampt_ff1[i, j], params.ksv[i, j], params.faif[i, j],
                    params.slope["gammaa"][i, j], params.infilt_limit[i, j], dt,
                )

        drained = 0.0
        for i in range(g.ny):
            for j in range(g.nx):
                if g.sink[i, j]:
                    drained += hs4[i, j] * g.area
                    hs4[i, j] = 0.0
        hr3 = dict(hr2)
        for (i, j) in g.down_riv:
            if g.sink[i, j]:
                channel_area = params.river["width"][(i, j)] * g.len_riv_grid[i, j]
                drained += hr3[(i, j)] * channel_area
                hr3[(i, j)] = 0.0

        diag = RefDiagnostics(qr=qr_end, qrs=qrs, aevp=aevp, infilt_rate=infilt_rate, drained_volume=drained)
        return hs4, gampt_ff2, hr3, diag

    def simulate(self, rain_seq, params: RefParams, dt: float, pet_seq=None, outlet_key: Optional[tuple] = None):
        hs, gampt_ff, hr = self.initial_state()
        qr_outlet = []
        hs_series = []
        hr_series = []
        T = rain_seq.shape[0]
        for t in range(T):
            pet_t = pet_seq[t] if pet_seq is not None else None
            hs, gampt_ff, hr, diag = self.step(hs, gampt_ff, hr, params, rain_seq[t], pet_t, dt)
            if outlet_key is not None:
                qr_outlet.append(diag.qr[outlet_key])
            hs_series.append(hs.copy())
            hr_series.append(dict(hr))
        out = {"hs_final": hs, "gampt_ff_final": gampt_ff, "hr_final": hr,
               "hs": np.stack(hs_series), "hr": hr_series}
        if outlet_key is not None:
            out["qr_outlet"] = np.array(qr_outlet)
        return out
