"""Vertical processes applied as an operator-split update after routing:
Green-Ampt infiltration (RRI_Infilt.f90: `infilt`) and evapotranspiration
(RRI_Evp.f90: `evp`). Both are simple per-cell elementwise updates, fully
vectorized and differentiable.
"""

from __future__ import annotations

import torch

_EPS = 1e-12


def infiltration(
    hs: torch.Tensor,
    gampt_ff: torch.Tensor,
    ksv: torch.Tensor,
    faif: torch.Tensor,
    gammaa: torch.Tensor,
    infilt_limit: torch.Tensor,
    dt: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Green-Ampt infiltration. `infilt_limit < 0` means "no limit"."""
    gampt_ff_safe = torch.clamp(gampt_ff, min=0.01)
    f = ksv * (1.0 + faif * gammaa / gampt_ff_safe)
    f = torch.minimum(f, hs / dt)
    limited = (infilt_limit >= 0.0) & (gampt_ff >= infilt_limit)
    f = torch.where(limited, torch.zeros_like(f), f)

    gampt_ff_new = gampt_ff + f * dt
    hs_new = torch.clamp(hs - f * dt, min=0.0)
    return hs_new, gampt_ff_new, f


def evaporation(
    hs: torch.Tensor,
    gampt_ff: torch.Tensor,
    pet: torch.Tensor,
    dt: float,
    evp_switch: int = 1,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Potential-ET extraction. `evp_switch = 1` also draws from the
    infiltrated store `gampt_ff` once `hs` is exhausted; `evp_switch = 2`
    (or any other value) draws from `hs` only.
    """
    cap = (hs + gampt_ff) if evp_switch == 1 else hs
    aevp = torch.minimum(pet, cap / dt)
    hs_new = hs - aevp * dt

    deficit = torch.clamp(hs_new, max=0.0)
    if evp_switch == 1:
        gampt_ff_new = torch.clamp(gampt_ff + deficit, min=0.0)
    else:
        gampt_ff_new = gampt_ff
    hs_new = torch.clamp(hs_new, min=0.0)
    return hs_new, gampt_ff_new, aevp
