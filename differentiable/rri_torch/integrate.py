"""Fixed-step explicit time integration.

The reference Fortran integrates both the slope and river ODEs with an
embedded, adaptive-step Runge-Kutta-Fehlberg 4(5) scheme (Cash-Karp
coefficients), independently sub-cycling within each output step until the
local error estimate is small enough (see the `b21..c6` coefficients and
`errmax`/`safety`/`pshrnk` logic in RRI.f90).

For a differentiable model built for gradient-based calibration, an
adaptive step count that changes with the (learnable) parameters is
undesirable: it makes the computational graph's shape parameter-dependent
and can introduce discontinuities in the loss surface at step-accept/
step-reject boundaries. We instead use a plain fixed-step RK4 integrator
with a user-chosen number of substeps per forcing interval -- accuracy is
controlled by picking `n_substeps` large enough for the fluxes/depths in
your domain (halve the step and re-run to check convergence), the same way
you would size `dt`/`dt_riv` in the original model.
"""

from __future__ import annotations

from typing import Callable

import torch

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


def integrate_fixed(rhs: RHS, y0: torch.Tensor, dt: float, n_substeps: int, clamp_min: float = 0.0) -> torch.Tensor:
    """Integrate dy/dt = rhs(y) over a total interval `dt`, split into
    `n_substeps` equal RK4 steps."""
    sub_dt = dt / n_substeps
    y = y0
    for _ in range(n_substeps):
        y = rk4_step(rhs, y, sub_dt, clamp_min=clamp_min)
    return y
