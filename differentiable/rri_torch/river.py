"""1-D river-network diffusive-wave routing.

Vectorized re-implementation of `funcr` / `qr_calc` / `hq_riv` in
RRI_Riv.f90, restricted to: a single D8 tree per catchment (no dams,
diversions, or discharge/water-level boundary files), and a default
rectangular channel (no `sec_map` cross-section tables).

River cells are stored as flat (n_riv,) tensors; the tree structure lives
in `Grid.down_riv` (index of the downstream river cell, or -1 at an
outlet) and is built once, non-differentiably, in `geometry.py`.
"""

from __future__ import annotations

from dataclasses import dataclass

import torch

from .geometry import Grid
from .hydraulics import hq_river

_EPS = 1e-12


@dataclass
class RiverParams:
    width: torch.Tensor       # (n_riv,) channel width [m]
    depth: torch.Tensor       # (n_riv,) channel depth below ground [m]
    ns_river: torch.Tensor    # scalar or (n_riv,) Manning's roughness [s/m^(1/3)]
    height: torch.Tensor      # (n_riv,) levee height above ground [m], 0 = no levee
    outlet_slope: float = 1e-3  # free-flow bed-slope assumed beyond the last outlet cell


def zb_riv(grid: Grid, params: RiverParams) -> torch.Tensor:
    """Riverbed elevation at each river cell: ground level minus channel depth."""
    return grid.gather_to_riv(grid.zb) - params.depth


def river_rhs(
    hr: torch.Tensor,
    grid: Grid,
    params: RiverParams,
) -> tuple[torch.Tensor, torch.Tensor]:
    """d(hr)/dt [m/s] for every river cell, and the discharge leaving each
    cell toward its downstream neighbour [m^3/s] (for diagnostics / the
    river-slope exchange step), following `funcr` in RRI_Riv.f90.
    """
    n = grid.n_riv
    if n == 0:
        return hr, hr

    down = grid.down_riv
    has_down = down >= 0
    safe_down = down.clamp(min=0)

    zb_p = zb_riv(grid, params)
    zb_n = torch.where(has_down, zb_p[safe_down], zb_p - grid.dis_riv * params.outlet_slope)

    hr_p = hr
    hr_n = torch.where(has_down, hr[safe_down], torch.zeros_like(hr))

    dh = ((zb_p + hr_p) - (zb_n + hr_n)) / grid.dis_riv

    hw_fwd = torch.where(zb_p >= zb_n, hr_p, torch.clamp(zb_p + hr_p - zb_n, min=0.0))
    hw_rev = torch.where(zb_n >= zb_p, hr_n, torch.clamp(zb_n + hr_n - zb_p, min=0.0))

    q_fwd = hq_river(hw_fwd, dh, params.width, params.ns_river)
    q_rev = hq_river(hw_rev, dh, params.width, params.ns_river)

    qr = torch.where(dh >= 0, q_fwd, -q_rev)

    inflow = torch.zeros_like(qr)
    if bool(has_down.any()):
        inflow = inflow.index_add(0, down[has_down], qr[has_down])

    channel_area = torch.clamp(params.width * grid.len_riv, min=_EPS)
    dhr_dt = (-qr + inflow) / channel_area
    return dhr_dt, qr
