"""Pre-fix ("legacy") versions of the two formulas fixed in HANDOFF.md
section 7.2, kept ONLY for the paper ablation in `ablation_solo.py` --
not used by the production model. Do not import this from `rri_torch`
itself; it exists purely to reconstruct "what the model computed before
each fix" for a controlled 2x2 comparison, via monkeypatching the
consuming modules' namespaces at ablation-run time.
"""
from __future__ import annotations

import torch

from rri_torch.tensor_ops import safe_sqrt

_G = 9.81
_MU1 = (2.0 / 3.0) ** 1.5
_MU2 = 0.35
_MU3 = 0.91
_EPS = 1e-9
_EPS_H = 1e-12


def legacy_hq_river(h: torch.Tensor, dh: torch.Tensor, width: torch.Tensor, ns_river: torch.Tensor) -> torch.Tensor:
    """Pre-fix `hq_river`: wide-channel approximation (R ~= h), i.e.
    `q = a * h**(5/3) * width` with no hydraulic-radius correction."""
    dh = dh.abs()
    a = safe_sqrt(dh) / torch.clamp(ns_river, min=_EPS_H)
    m = 5.0 / 3.0
    return a * torch.clamp(h, min=0.0) ** m * width


def _legacy_weir_flow(h1, h2, mu2, mu3, dt, len_riv, area):
    """Pre-fix `_weir_flow`: missing the "both banks" factor of 2."""
    free = mu2 * h1 * safe_sqrt(2.0 * _G * h1) * dt * len_riv / area
    submerged = mu3 * h2 * safe_sqrt(2.0 * _G * (h1 - h2)) * dt * len_riv / area
    q = torch.where(torch.clamp(h2, min=0.0) / torch.clamp(h1, min=_EPS) <= 2.0 / 3.0, free, submerged)
    return torch.where(h1 <= 0.0, torch.zeros_like(h1), q)


def legacy_river_slope_exchange(hr, hs, geom, dt):
    """Pre-fix `river_slope_exchange`: identical to the current
    `exchange.river_slope_exchange` except `hrs_a` and `_weir_flow`'s
    `free`/`submerged` terms are missing the `* 2.0` both-banks factor
    (see `exchange._weir_flow`'s docstring for what that factor is)."""
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

    hrs_a = torch.clamp(_MU1 * hs_top * safe_sqrt(_G * hs_top) * dt * len_riv / area, max=hs_top)

    h1_c = hr_top - height
    h2_c = hs_top - height
    hrs_c = -_legacy_weir_flow(h1_c, h2_c, _MU2, _MU3, dt, len_riv, area)
    cap_c = -(hr_top - height) * (len_riv * width / area)
    hrs_c = torch.where(hrs_c.abs() > cap_c.abs(), cap_c, hrs_c)

    h1_d = hs_top - height
    h2_d = hr_top - height
    hrs_d = _legacy_weir_flow(h1_d, h2_d, _MU2, _MU3, dt, len_riv, area)
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

    needs_corr_ad = (case_a | case_d) & (hr_new - geom.depth >= -1e-5) & (hr_new - geom.depth > hs_new)
    needs_corr_c = case_c & (hs_new > hr_new - geom.depth)
    from rri_torch.exchange import N_CORRECTION
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
