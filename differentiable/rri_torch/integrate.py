"""Explicit time integration: a fixed-step RK4 integrator, and a real
adaptive-step Runge-Kutta-Fehlberg 4(5) integrator (Cash-Karp
coefficients) matching RRI.f90's own scheme exactly.

The reference Fortran integrates both the slope and river ODEs with an
embedded, adaptive-step RKF45 scheme, independently sub-cycling within
each output step until the local error estimate is small enough (see the
`b21..c6` coefficients and `errmax`/`safety`/`pshrnk` logic in RRI.f90,
and `reference/integrate_ref.py`'s already-verified NumPy port of it).

For a differentiable model built for gradient-based calibration, an
adaptive step count that changes with the (learnable) parameters is
undesirable in general: it makes the computational graph's shape
parameter-dependent, and can introduce discontinuities in the loss
surface at step-accept/step-reject boundaries. `integrate_fixed` sidesteps
this entirely with a user-chosen, parameter-independent substep count.

`integrate_adaptive` below trades that guarantee for real accuracy control
(useful for validation against the reference Fortran, e.g. on real-world
scenarios whose stiffness varies a lot in space -- see
`differentiable/HANDOFF.md` section 6a for a case where a fixed substep
count needed manual tuning per reach). It IS still usable under autograd
(gradients flow through the arithmetic of every accepted RK stage, exactly
like `torchdiffeq`'s default direct-backprop mode) -- the only
non-differentiable part is the accept/reject and step-size decision
itself, which reads the error estimate through `.item()`. This means a
gradient step can change *which* steps get taken without that change
itself producing a gradient signal, so the loss surface can have small
kinks at step boundaries; for calibration work prefer `integrate_fixed`
with a substep count picked to resolve the domain's stiffness, and reserve
`integrate_adaptive` for validation/forward-only runs where exact
agreement with the reference Fortran matters more than a perfectly smooth
loss surface.
"""

from __future__ import annotations

from typing import Callable, NamedTuple, Optional

import torch
import torch.utils.checkpoint

RHS = Callable[[torch.Tensor], torch.Tensor]


def rk4_step(rhs: RHS, y0: torch.Tensor, dt: float, clamp_min: float = 0.0) -> torch.Tensor:
    """One classical RK4 step for dy/dt = rhs(y), with intermediate and
    final states clamped to `clamp_min` (RRI clamps negative depths to
    zero after every Runge-Kutta stage)."""
    def c(y):
        return torch.clamp(y, min=clamp_min)

    k1 = rhs(c(y0))
    k2 = rhs(c(y0 + 0.5 * dt * k1))
    k3 = rhs(c(y0 + 0.5 * dt * k2))
    k4 = rhs(c(y0 + dt * k3))
    y1 = y0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return c(y1)


def integrate_fixed(
    rhs: RHS, y0: torch.Tensor, dt: float, n_substeps: int, clamp_min: float = 0.0,
    on_step: Optional[Callable[[torch.Tensor, float], None]] = None,
    checkpoint_substeps: bool = False,
) -> torch.Tensor:
    """Integrate dy/dt = rhs(y) over a total interval `dt`, split into
    `n_substeps` equal RK4 steps. If given, `on_step(y_new, sub_dt)` is
    called after every substep -- e.g. to accumulate a `ddt`-weighted time
    average of some quantity derived from the state, analogous to RRI's
    own `qr_ave` (see `RRI.f90` lines 644-715 and HANDOFF.md section 6a).

    `checkpoint_substeps=True` wraps each individual RK4 substep in its
    own `torch.utils.checkpoint.checkpoint`, on top of whatever
    step-level checkpointing the caller may already be using (nested
    checkpointing composes correctly under `use_reentrant=False`). This
    matters because a single call already unrolls `n_substeps` RK4 steps
    (4 RHS evals each) into one autograd graph -- on a real (not toy)
    grid with `n_substeps` in the tens to hundreds (needed for fixed-step
    stability, see HANDOFF.md section 6a), *that one call's* local
    backward graph is itself too large to hold in GPU memory even when
    `RRIModel.simulate(use_checkpointing=True)` already checkpoints at
    the outer time-step level (confirmed: OOMs at 25+ GB on the real
    Solo scenario's 68544-cell grid at n_substeps_slope=40, *independent
    of how many outer steps are being backpropagated* -- see
    `test_checkpoint_substeps_match_uncheckpointed`). Checkpointing each
    substep individually bounds the resident graph to one RK4 step (4
    RHS evals) at a time instead of all `n_substeps` of them.

    Incompatible with `on_step`: `on_step` has side effects (it mutates
    an accumulator via a closure, e.g. `RRIModel.step`'s `qr_avg`
    tracking), and checkpoint recomputation would silently re-run those
    side effects a second time during backward if `on_step` sat inside
    the checkpointed region. Since the only current `on_step` use
    (`track_qr_avg`) is a validation-only diagnostic that calibration
    runs (the only place `checkpoint_substeps=True` is needed) don't use,
    this is simply disallowed rather than made "correct but surprising".
    """
    if checkpoint_substeps and on_step is not None:
        raise ValueError("checkpoint_substeps=True is incompatible with on_step (see docstring)")
    sub_dt = dt / n_substeps
    y = y0
    for _ in range(n_substeps):
        if checkpoint_substeps:
            y = torch.utils.checkpoint.checkpoint(
                rk4_step, rhs, y, sub_dt, clamp_min, use_reentrant=False,
            )
        else:
            y = rk4_step(rhs, y, sub_dt, clamp_min=clamp_min)
        if on_step is not None:
            on_step(y, sub_dt)
    return y


# ---------------------------------------------------------------------
# Adaptive Cash-Karp RKF45, matching RRI_Mod2.f90 / RRI.f90 exactly.
# Coefficients transcribed from reference/integrate_ref.py (already
# cross-checked against the Fortran's own `runge_mod` in
# tests/test_reference_agreement.py).
# ---------------------------------------------------------------------

_A2, _A3, _A4, _A5, _A6 = 0.2, 0.3, 0.6, 1.0, 0.875
_B21 = 0.2
_B31, _B32 = 3.0 / 40.0, 9.0 / 40.0
_B41, _B42, _B43 = 0.3, -0.9, 1.2
_B51, _B52, _B53, _B54 = -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0
_B61, _B62, _B63, _B64, _B65 = (
    1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0,
)
_C1, _C3, _C4, _C6 = 37.0 / 378.0, 250.0 / 621.0, 125.0 / 594.0, 512.0 / 1771.0
_DC1 = _C1 - 2825.0 / 27648.0
_DC3 = _C3 - 18575.0 / 48384.0
_DC4 = _C4 - 13525.0 / 55296.0
_DC5 = -277.0 / 14336.0
_DC6 = _C6 - 0.25

_SAFETY, _PSHRNK = 0.9, -0.25


class AdaptiveStats(NamedTuple):
    n_accepted: int
    n_rejected: int
    min_ddt: float
    max_ddt: float


def integrate_adaptive(
    rhs: RHS,
    y0: torch.Tensor,
    dt_outer: float,
    eps: float,
    ddt_min: float,
    mask: Optional[torch.Tensor] = None,
    clamp_min: float = 0.0,
    max_iters: int = 100_000,
    on_step: Optional[Callable[[torch.Tensor, float], None]] = None,
) -> tuple[torch.Tensor, AdaptiveStats]:
    """Integrate dy/dt = rhs(y) over `dt_outer` with RRI's real adaptive
    Cash-Karp RKF45 step control: shrink a single shared step size `ddt`
    whenever the embedded 4th/5th-order error estimate (relative to `eps`,
    RRI's absolute-storage tolerance -- see HANDOFF.md on why this needs
    to be sized to your domain's typical depths) exceeds 1, accept
    otherwise, floor `ddt` at `ddt_min` (matching RRI's own "give up
    shrinking, accept anyway" floor).

    `mask`, if given, restricts the error estimate to those entries (e.g.
    a domain mask for slope cells) -- inactive entries can't veto the
    step size. All arithmetic is ordinary differentiable tensor ops;
    only the step-size/accept-reject decision reads `errmax` through
    `.item()` (see module docstring for the autograd implications).

    `on_step(y_new, ddt)`, if given, is called after every *accepted*
    substep (not rejected ones) -- see `integrate_fixed`'s docstring for
    why you'd want this (a `qr_ave`-style time average).
    """
    def c(y):
        return torch.clamp(y, min=clamp_min)

    time = 0.0
    y = y0
    ddt = dt_outer
    n_accepted = 0
    n_rejected = 0
    min_ddt = dt_outer
    max_ddt = 0.0

    for _ in range(max_iters):
        if time + ddt > dt_outer:
            ddt = dt_outer - time

        while True:
            k1 = rhs(c(y))
            k2 = rhs(c(y + _B21 * ddt * k1))
            k3 = rhs(c(y + ddt * (_B31 * k1 + _B32 * k2)))
            k4 = rhs(c(y + ddt * (_B41 * k1 + _B42 * k2 + _B43 * k3)))
            k5 = rhs(c(y + ddt * (_B51 * k1 + _B52 * k2 + _B53 * k3 + _B54 * k4)))
            k6 = rhs(c(y + ddt * (_B61 * k1 + _B62 * k2 + _B63 * k3 + _B64 * k4 + _B65 * k5)))

            y_new = c(y + ddt * (_C1 * k1 + _C3 * k3 + _C4 * k4 + _C6 * k6))
            err = ddt * (_DC1 * k1 + _DC3 * k3 + _DC4 * k4 + _DC5 * k5 + _DC6 * k6)
            if mask is not None:
                err = torch.where(mask, err, torch.zeros_like(err))
            errmax = (err.abs().max() / eps).item()

            if errmax > 1.0 and ddt > ddt_min:
                n_rejected += 1
                ddt = max(_SAFETY * ddt * errmax ** _PSHRNK, 0.5 * ddt)
                continue
            else:
                if time + ddt > dt_outer:
                    ddt = dt_outer - time
                time += ddt
                y = y_new
                n_accepted += 1
                min_ddt = min(min_ddt, ddt)
                max_ddt = max(max_ddt, ddt)
                if on_step is not None:
                    on_step(y, ddt)
                break

        if time >= dt_outer:
            break
    else:
        raise RuntimeError(
            f"integrate_adaptive did not reach dt_outer={dt_outer} within {max_iters} outer iterations "
            f"(reached time={time}); the domain may be too stiff for this eps/ddt_min combination."
        )

    return y, AdaptiveStats(n_accepted=n_accepted, n_rejected=n_rejected, min_ddt=min_ddt, max_ddt=max_ddt)
