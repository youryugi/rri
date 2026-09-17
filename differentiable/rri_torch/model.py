"""Assembles slope routing, river routing, river-slope exchange,
infiltration and evaporation into one time-stepping model, preserving the
operator-splitting order used by the main loop in RRI.f90:

  1. river routing over the full step dt (RK, ignoring slope exchange)
  2. slope routing over the full step dt (RK, rainfall forcing only)
  3. evapotranspiration (explicit, post-routing)
  4. river <-> slope exchange (weir formulas)
  5. infiltration (Green-Ampt, explicit, post-exchange)
  6. domain-boundary ("sink") cells drain to zero

Not implemented (see README for the full list and rationale): groundwater,
dams, diversions, levee breaks, sediment transport, boundary-condition
files, and non-rectangular cross-sections. Adaptive RKF45 time-stepping
*is* available (see `integrate.integrate_adaptive` and the `adaptive=`
flag below) but defaults to off -- fixed-step RK4 remains the default
because it keeps the computational graph's shape parameter-independent,
which matters for gradient-based calibration (see `integrate.py`'s
module docstring for the full tradeoff).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import torch
import torch.utils.checkpoint

from .geometry import Grid
from .slope import SlopeParams, slope_rhs
from .river import RiverParams, river_rhs
from .exchange import ExchangeGeometry, river_slope_exchange
from .processes import infiltration, evaporation
from .integrate import integrate_fixed, integrate_adaptive, AdaptiveStats


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
    river_stats: Optional[AdaptiveStats] = None   # only set when adaptive=True
    slope_stats: Optional[AdaptiveStats] = None   # only set when adaptive=True
    qr_avg: Optional[torch.Tensor] = None         # (n_riv,), only set when track_qr_avg=True


class RRIModel:
    def __init__(
        self,
        grid: Grid,
        n_substeps_slope: int = 4,
        n_substeps_river: int = 4,
        adaptive: bool = False,
        eps: float = 0.01,
        ddt_min_slope: float = 1.0,
        ddt_min_river: float = 0.1,
        track_qr_avg: bool = False,
        checkpoint_substeps: bool = False,
    ):
        """`adaptive=True` switches river and slope routing to RRI's real
        adaptive Cash-Karp RKF45 stepping (`integrate.integrate_adaptive`)
        instead of fixed-substep RK4; `n_substeps_slope`/`n_substeps_river`
        are then ignored. `eps`/`ddt_min_slope`/`ddt_min_river` default to
        RRI's own hardcoded constants (`RRI_Mod2.f90`: eps=0.01,
        ddt_min_slo=1.0, ddt_min_riv=0.1) -- eps is an *absolute* storage
        tolerance [m], not relative, so it needs resizing for domains with
        very different typical depths than a real river network's (see
        HANDOFF.md's note on this for a small synthetic catchment).

        `track_qr_avg=True` additionally reports `StepDiagnostics.qr_avg`:
        the river discharge *time-averaged* over the outer step (a
        ddt-weighted Riemann sum over every accepted substep's end state),
        matching what RRI's own `hydro.txt` output actually is (`qr_ave`
        in `RRI.f90`, accumulated as `qr_idx * ddt` and normalized by
        `dt` -- see HANDOFF.md section 6a) as opposed to the instantaneous
        end-of-step value `StepDiagnostics.qr` always reports. Off by
        default since it costs one extra `river_rhs` call per substep.

        `checkpoint_substeps=True` (only meaningful with `adaptive=False`)
        additionally checkpoints every individual RK4 substep inside
        `step`'s river/slope `integrate_fixed` calls -- needed alongside
        `simulate(use_checkpointing=True)` on real (not toy) grids, where
        `n_substeps_{slope,river}` in the tens to hundreds (see
        HANDOFF.md section 6a) makes even *one* outer step's local
        backward graph too large to hold, independent of how many outer
        steps are being backpropagated (see `integrate_fixed`'s
        docstring). Incompatible with `track_qr_avg=True` (see
        `integrate_fixed`'s `on_step` note) -- checked eagerly here."""
        if checkpoint_substeps and track_qr_avg:
            raise ValueError("checkpoint_substeps=True is incompatible with track_qr_avg=True")
        self.grid = grid
        self.n_substeps_slope = n_substeps_slope
        self.n_substeps_river = n_substeps_river
        self.adaptive = adaptive
        self.eps = eps
        self.ddt_min_slope = ddt_min_slope
        self.ddt_min_river = ddt_min_river
        self.track_qr_avg = track_qr_avg
        self.checkpoint_substeps = checkpoint_substeps

    def initial_state(self, dtype: torch.dtype = torch.float64):
        g = self.grid
        device = g.domain.device
        hs = torch.zeros((g.ny, g.nx), dtype=dtype, device=device)
        gampt_ff = torch.zeros((g.ny, g.nx), dtype=dtype, device=device)
        hr = torch.zeros((g.n_riv,), dtype=dtype, device=device)
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
        river_stats = slope_stats = None

        # 1. river routing (independent of slope for the duration of dt)
        #
        # `qr_stage_log`, when track_qr_avg is on, records the qr byproduct
        # of every internal RHS evaluation (i.e. every RK stage: 4 for
        # fixed RK4, 6 for adaptive Cash-Karp) -- exactly what RRI.f90
        # accumulates into `qr_ave_idx` (lines 644-715: sum the 6 stage
        # qr's, divide by 6, weight by ddt, sum over substeps, divide by
        # dt). A rejected-then-retried substep re-evaluates all stages
        # from scratch, so `on_river_step` only needs the *last*
        # `n_stages` log entries -- whatever an earlier rejected attempt
        # logged is simply superseded, never included.
        n_stages = 6 if self.adaptive else 4
        qr_stage_log: list[torch.Tensor] = []

        def river_rhs_fn(hr_):
            dhrdt, qr_stage = river_rhs(hr_, g, params.river)
            if self.track_qr_avg:
                qr_stage_log.append(qr_stage)
            return dhrdt

        qr_avg = None
        on_river_step = None
        if self.track_qr_avg:
            qr_accum = torch.zeros(g.n_riv, dtype=hr.dtype, device=hr.device)
            ddt_accum = 0.0

            def on_river_step(hr_new, ddt):
                nonlocal qr_accum, ddt_accum
                stages = qr_stage_log[-n_stages:]
                qr_accum = qr_accum + (sum(stages) / len(stages)) * ddt
                ddt_accum += ddt
                qr_stage_log.clear()

        if self.adaptive:
            hr1, river_stats = integrate_adaptive(
                river_rhs_fn, hr, dt, self.eps, self.ddt_min_river, on_step=on_river_step,
            )
        else:
            hr1 = integrate_fixed(
                river_rhs_fn, hr, dt, self.n_substeps_river, on_step=on_river_step,
                checkpoint_substeps=self.checkpoint_substeps,
            )
        _, qr_end = river_rhs(hr1, g, params.river)
        if self.track_qr_avg:
            qr_avg = qr_accum / ddt_accum

        # 2. slope routing (rainfall forcing only)
        def slope_rhs_fn(hs_):
            dhsdt, _ = slope_rhs(hs_, g.zb, g.domain, params.slope, rain, g.dx, g.dy, g.area)
            return dhsdt

        if self.adaptive:
            hs1, slope_stats = integrate_adaptive(
                slope_rhs_fn, hs, dt, self.eps, self.ddt_min_slope, mask=g.domain,
            )
        else:
            hs1 = integrate_fixed(
                slope_rhs_fn, hs, dt, self.n_substeps_slope,
                checkpoint_substeps=self.checkpoint_substeps,
            )

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

        diag = StepDiagnostics(qr=qr_end, qrs=qrs, aevp=aevp, infilt_rate=infilt_rate, drained_volume=drained,
                                river_stats=river_stats, slope_stats=slope_stats, qr_avg=qr_avg)
        return hs4, gampt_ff2, hr3, diag

    def _step_flat(self, hs, gampt_ff, hr, params, rain_t, pet_t, dt):
        """Same as `step`, but returns a flat tuple of tensors only (no
        `StepDiagnostics` dataclass) -- `torch.utils.checkpoint.checkpoint`
        needs its wrapped function's output to be tensors (or nested
        tensor containers) it can attach backward hooks to, not an
        arbitrary dataclass. Used by `simulate(..., use_checkpointing=True)`.
        """
        hs4, gampt_ff2, hr3, diag = self.step(hs, gampt_ff, hr, params, rain_t, pet_t, dt)
        qr_avg = diag.qr_avg if diag.qr_avg is not None else torch.zeros(0, dtype=hr.dtype, device=hr.device)
        return hs4, gampt_ff2, hr3, diag.qr, diag.drained_volume, qr_avg

    def simulate(
        self,
        rain_seq: torch.Tensor,
        params: RRIParams,
        dt: float,
        pet_seq: Optional[torch.Tensor] = None,
        outlet_riv_index: Optional[int] = None,
        init_state: Optional[tuple] = None,
        record_full_state: bool = False,
        use_checkpointing: bool = False,
    ) -> dict:
        """Run the model over a forcing time series.

        rain_seq: (T, ny, nx) rainfall rate [m/s] for each step.
        pet_seq: optional (T, ny, nx) potential ET rate [m/s].
        outlet_riv_index: river-cell index to report a discharge series
            for (typically the catchment outlet); if None, no `qr_outlet`
            key is returned.
        record_full_state: if True, also stack the full `hs`/`hr` history
            (memory-heavy for long runs / big grids).
        use_checkpointing: if True, wrap each outer step in
            `torch.utils.checkpoint.checkpoint` -- trades ~2x forward
            compute (each step's internal substeps get recomputed once
            during the backward pass) for O(n_substeps) -> O(1) memory per
            step, since none of a step's *internal* substep activations
            need to stay resident between steps, only the state handed to
            the next one. Backpropagating through more than a couple
            hundred un-checkpointed fixed-RK4 substeps on a real (not
            toy) grid exhausts GPU memory outright (confirmed: a 32GB GPU
            OOMs backpropagating just 200 outer steps at
            n_substeps_slope=40 on the 18582-cell Solo grid) -- turn this
            on for any calibration run longer than a handful of steps.
            Forward-only (`torch.no_grad()`) runs get no benefit from
            this and should leave it off.
        """
        hs, gampt_ff, hr = init_state if init_state is not None else self.initial_state(dtype=rain_seq.dtype)

        qr_outlet_series = []
        qr_avg_outlet_series = []
        drained_series = []
        hs_series = []
        hr_series = []

        T = rain_seq.shape[0]
        for t in range(T):
            pet_t = pet_seq[t] if pet_seq is not None else None
            if use_checkpointing:
                hs, gampt_ff, hr, qr, drained, qr_avg = torch.utils.checkpoint.checkpoint(
                    self._step_flat, hs, gampt_ff, hr, params, rain_seq[t], pet_t, dt, use_reentrant=False,
                )
                qr_avg = qr_avg if qr_avg.numel() > 0 else None
            else:
                hs, gampt_ff, hr, diag = self.step(hs, gampt_ff, hr, params, rain_seq[t], pet_t, dt)
                qr, drained, qr_avg = diag.qr, diag.drained_volume, diag.qr_avg
            drained_series.append(drained)
            if outlet_riv_index is not None:
                qr_outlet_series.append(qr[outlet_riv_index])
                if qr_avg is not None:
                    qr_avg_outlet_series.append(qr_avg[outlet_riv_index])
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
            if qr_avg_outlet_series:
                out["qr_avg_outlet"] = torch.stack(qr_avg_outlet_series)
        if record_full_state:
            out["hs"] = torch.stack(hs_series)
            out["hr"] = torch.stack(hr_series)
        return out
