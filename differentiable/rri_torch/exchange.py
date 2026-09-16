"""River <-> slope water exchange (weir-type formulas).

Vectorized re-implementation of `funcrs_dt` in RRI_RivSlo.f90, restricted
to the default 4-case weir exchange (no dam handling, no
`river_overtop_neighbor_switch` spreading to neighbouring cells, no
`sec_map` cross-section tables so channel top width is just `width`).

The original Fortran resolves the "don't overshoot equilibrium" correction
with a `do ... exit` loop of up to 10 iterations per cell. Since the
iteration count would otherwise depend on data (breaking the static graph
shape needed for autograd), we replace it with a fixed number of
correction passes; each pass is a no-op wherever the cell has already
converged, so this is exact once `N_CORRECTION` is large enough and a
close approximation otherwise (see README).
"""

from __future__ import annotations

from dataclasses import dataclass

import torch

from .tensor_ops import safe_sqrt

N_CORRECTION = 6
_G = 9.81
_MU1 = (2.0 / 3.0) ** 1.5
_MU2 = 0.35
_MU3 = 0.91
_EPS = 1e-9


@dataclass
class ExchangeGeometry:
    riv_mask: torch.Tensor   # bool (ny, nx)
    depth: torch.Tensor      # float (ny, nx), channel depth below ground (0 off-river)
    height: torch.Tensor     # float (ny, nx), levee height above ground (0 off-river)
    width: torch.Tensor      # float (ny, nx), channel top width (0 off-river)
    len_riv: torch.Tensor    # float (ny, nx), channel reach length (0 off-river)
    area: float


def _weir_flow(h1: torch.Tensor, h2: torch.Tensor, mu2: float, mu3: float, dt: float, len_riv: torch.Tensor, area: float) -> torch.Tensor:
    """Free/submerged broad-crested weir formula shared by cases (c)/(d)."""
    free = mu2 * h1 * safe_sqrt(2.0 * _G * h1) * dt * len_riv / area
    submerged = mu3 * h2 * safe_sqrt(2.0 * _G * (h1 - h2)) * dt * len_riv / area
    q = torch.where(torch.clamp(h2, min=0.0) / torch.clamp(h1, min=_EPS) <= 2.0 / 3.0, free, submerged)
    return torch.where(h1 <= 0.0, torch.zeros_like(h1), q)


def river_slope_exchange(
    hr: torch.Tensor,
    hs: torch.Tensor,
    geom: ExchangeGeometry,
    dt: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Exchange water between the river state `hr` and slope state `hs`
    over a time interval `dt`.

    Returns (hr_new, hs_new, qrs) where qrs [m/s] is the net rate from
    slope to river (positive = slope -> river), matching RRI's `qrs`
    output convention.
    """
    riv = geom.riv_mask
    height = torch.clamp(geom.height, min=0.0)
    width = geom.width
    len_riv = geom.len_riv
    area = geom.area

    hr_top = hr - geom.depth
    hs_top = hs

    case_a = riv & (
        ((height == 0.0) & (hr_top < 0.0))
        | ((height > 0.0) & (hr_top < 0.0) & (hs_top <= height))
    )
    case_b = riv & (height > 0.0) & (hs_top <= height) & (hr_top <= height) & (hr_top >= 0.0)
    case_c = riv & (hs_top <= hr_top) & (hr_top >= height) & ~case_a & ~case_b
    case_d = riv & (hs_top >= hr_top) & (hs_top >= height) & ~case_a & ~case_b & ~case_c
    # any remaining river cell (numerically ambiguous boundary between
    # cases) exchanges nothing rather than raising, unlike the Fortran.

    hrs_a = torch.clamp(_MU1 * hs_top * safe_sqrt(_G * hs_top) * dt * len_riv / area, max=hs_top)

    h1_c = hr_top - height
    h2_c = hs_top - height
    hrs_c = -_weir_flow(h1_c, h2_c, _MU2, _MU3, dt, len_riv, area)
    cap_c = -(hr_top - height) * (len_riv * width / area)
    hrs_c = torch.where(hrs_c.abs() > cap_c.abs(), cap_c, hrs_c)

    h1_d = hs_top - height
    h2_d = hr_top - height
    hrs_d = _weir_flow(h1_d, h2_d, _MU2, _MU3, dt, len_riv, area)
    hrs_d = torch.where(hrs_d > (hs_top - height), hs_top - height, hrs_d)

    hrs = torch.zeros_like(hs)
    hrs = torch.where(case_d, hrs_d, hrs)
    hrs = torch.where(case_c, hrs_c, hrs)
    hrs = torch.where(case_b, torch.zeros_like(hs), hrs)
    hrs = torch.where(case_a, hrs_a, hrs)
    hrs = hrs * riv

    channel_area = torch.clamp(width * len_riv, min=_EPS)

    def apply(hs_, hr_, hrs_):
        hs_ = hs_ - hrs_
        hr_ = hr_ + hrs_ * area / channel_area
        return hs_, hr_

    hs_new, hr_new = apply(hs, hr, hrs)

    # Fixed-count correction passes to prevent hr/hs from crossing past
    # equilibrium, replacing the Fortran `do count = 1, 10 ... exit` loop.
    needs_corr_ad = (case_a | case_d) & (hr_new - geom.depth >= -1e-5) & (hr_new - geom.depth > hs_new)
    needs_corr_c = case_c & (hs_new > hr_new - geom.depth)
    for _ in range(N_CORRECTION):
        hs_top_i = hs_new
        hr_top_i = hr_new - geom.depth
        mask = (needs_corr_ad | needs_corr_c) & riv
        ar = len_riv * width / area
        corr = (hs_top_i - hr_top_i) / (1.0 + 1.0 / torch.clamp(ar, min=_EPS))
        corr = corr * mask
        hs_new = hs_new - corr
        hr_new = hr_new + corr * area / channel_area
        hrs = hrs + corr
        needs_corr_ad = (case_a | case_d) & (hr_new - geom.depth >= -1e-5) & (hr_new - geom.depth > hs_new) & ((hs_top_i - hr_top_i).abs() > 1e-6)
        needs_corr_c = case_c & (hs_new > hr_new - geom.depth) & ((hs_top_i - hr_top_i).abs() > 1e-6)

    hs_new = torch.where(riv, hs_new, hs)
    hr_new = torch.where(riv, hr_new, hr)
    qrs = torch.where(riv, hrs / dt, torch.zeros_like(hrs))
    return hr_new, hs_new, qrs
