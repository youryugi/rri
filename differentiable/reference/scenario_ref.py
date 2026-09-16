"""Builds the reference (RefGrid, RefParams) from the exact same raw
scenario dict used by `examples.synthetic_basin.build_synthetic_basin`,
so the differentiable model and this independent loop-based port run on
provably identical catchments/parameters.
"""

from __future__ import annotations

import numpy as np

from examples.synthetic_basin import build_raw_scenario
from .geometry_ref import build_ref_grid
from .model_ref import RefParams


def build_synthetic_basin_ref(**kwargs):
    """Returns (grid, params, outlet_key) for the reference port."""
    raw = build_raw_scenario(**kwargs)
    ny, nx = raw["ny"], raw["nx"]

    grid = build_ref_grid(
        domain=raw["domain"], zb=raw["zb"], dir_grid=raw["dir_grid"],
        riv_mask=raw["riv_mask"], len_riv=raw["len_riv"],
        dx=raw["dx"], dy=raw["dy"],
    )

    def full(value):
        return np.full((ny, nx), value, dtype=np.float64)

    slope = dict(
        ns=full(raw["ns_slope"]), ka=full(raw["ka"]), beta=full(raw["beta"]),
        da=full(raw["da"]), dm=full(raw["dm"]),
        soildepth=full(raw["soildepth"]), gammaa=full(raw["gammaa"]),
        min_wc4latflow=full(raw["min_wc4latflow"]),
    )

    river = dict(
        width={k: raw["river_width"] for k in grid.down_riv},
        depth={k: raw["river_depth"] for k in grid.down_riv},
        height={k: raw["river_height"] for k in grid.down_riv},
        ns_river=raw["ns_river"],
        outlet_slope=raw["outlet_slope"],
    )

    params = RefParams(
        slope=slope, river=river,
        ksv=full(raw["ksv"]), faif=full(raw["faif"]), infilt_limit=full(raw["infilt_limit"]),
        evp_switch=0,
    )

    outlet_key = max(grid.down_riv.keys(), key=lambda ij: ij[0])  # southernmost river cell
    return grid, params, outlet_key
