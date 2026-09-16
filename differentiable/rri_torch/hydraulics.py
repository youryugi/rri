"""Vectorized, differentiable versions of the storage <-> discharge and
water-level relationships in RRI_Slope.f90 (`hq`, `h2lev`) and
RRI_Riv.f90 (`hq_riv`), restricted to the default rectangular channel
(no `sec_map` cross-section lookup tables -- see README for scope notes).

All functions are elementwise over arbitrarily-shaped tensors and are safe
to call with `requires_grad=True` parameters: branch selection uses
`torch.where` so both branches are evaluated, and denominators that could
be exactly zero on the "wrong" branch are clamped before division so no
NaN gradients leak through the unused branch.
"""

from __future__ import annotations

import torch

from .tensor_ops import safe_sqrt

_EPS = 1e-12


def h2lev(h: torch.Tensor, soildepth: torch.Tensor, gammaa: torch.Tensor) -> torch.Tensor:
    """Storage depth `h` -> actual water level `lev` measured from the
    slope-cell's ground surface datum used for the diffusive-wave slope
    term (RRI_Slope.f90: `h2lev`).

    When h is below the soil's saturated-storage capacity (soildepth *
    gammaa), the water is sitting unsaturated within the soil profile and
    `lev` is inflated by 1/gammaa (porosity) to get the true water-table
    position. Above that capacity, extra storage sits as 1:1 surface
    ponding on top of the soil column.
    """
    da_temp = soildepth * gammaa
    safe_gammaa = torch.clamp(gammaa, min=_EPS)
    lev_unsat = h / safe_gammaa
    lev_sat = soildepth + (h - da_temp)
    return torch.where(h >= da_temp, lev_sat, lev_unsat)


def hq_slope(
    ns: torch.Tensor,
    ka: torch.Tensor,
    da: torch.Tensor,
    dm: torch.Tensor,
    beta: torch.Tensor,
    h: torch.Tensor,
    dh: torch.Tensor,
    edge_len: torch.Tensor,
    area: float,
    soildepth: torch.Tensor,
    gammaa: torch.Tensor,
    min_wc4latflow: torch.Tensor,
) -> torch.Tensor:
    """Lateral discharge (per unit slope-cell area, i.e. an equivalent
    depth rate [m/s]) driven from a cell with storage depth `h` and
    water-surface gradient magnitude `dh` (already >= 0), following RRI's
    tri-linear two-layer storage-discharge law (RRI_Slope.f90: `hq`):

      h <  dm : nonlinear unsaturated subsurface flow, exponent `beta`
      dm<=h<da: linear (Darcy) saturated subsurface flow
      h >= da : saturated subsurface flow + Manning overland flow

    `edge_len` is the cross-flow cell-edge length used to convert the
    per-unit-width discharge into a flux normalized by cell area.

    Deviation from the reference Fortran: the original `hq()` only
    suppresses flow below `min_hs` (a minimum wetness before lateral
    subsurface flow can start) when it is called for the *downhill*
    (outflow) side of an edge; the uphill/inflow-side call is passed the
    signed (negative) gradient and the cutoff can never trigger there. We
    apply the cutoff symmetrically to whichever cell is providing water
    (see `slope.py`), which is the physically consistent choice for this
    reimplementation.
    """
    dh = dh.abs()
    min_hs = torch.where(ka > 0, soildepth * gammaa * min_wc4latflow, torch.zeros_like(ka))

    km = torch.where(beta > 0, ka / torch.clamp(beta, min=_EPS), torch.zeros_like(ka))
    vm = km * dh
    va = torch.where(da > 0, ka * dh, torch.zeros_like(ka))
    al = safe_sqrt(dh) / torch.clamp(ns, min=_EPS)
    m = 5.0 / 3.0

    dm_safe = torch.clamp(dm, min=_EPS)
    da_safe = torch.clamp(da, min=_EPS)
    h_safe = torch.clamp(h, min=0.0)

    ratio_safe = torch.clamp(h_safe / dm_safe, min=_EPS)
    q_unsat = vm * dm * ratio_safe ** beta
    q_sat = vm * dm + va * (h_safe - dm)
    q_surf = vm * dm + va * (da_safe - dm) + al * torch.clamp(h_safe - da_safe, min=0.0) ** m

    q = torch.where(h_safe < dm, q_unsat, torch.where(h_safe < da, q_sat, q_surf))
    q = torch.where(h_safe <= min_hs, torch.zeros_like(q), q)
    return q * edge_len / area


def hq_river(h: torch.Tensor, dh: torch.Tensor, width: torch.Tensor, ns_river: torch.Tensor) -> torch.Tensor:
    """Manning's-equation discharge [m^3/s] for a wide rectangular channel
    (RRI_Riv.f90: `hq_riv`, without the optional cross-section table).
    """
    dh = dh.abs()
    a = safe_sqrt(dh) / torch.clamp(ns_river, min=_EPS)
    m = 5.0 / 3.0
    return a * torch.clamp(h, min=0.0) ** m * width
