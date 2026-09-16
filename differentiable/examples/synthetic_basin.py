"""Builds a small synthetic V-shaped catchment for testing and the
calibration demo: a straight river channel running down the middle
column, draining two symmetric hillslopes, exiting at the bottom row.

There is no real DEM/rainfall/discharge data bundled with this
repository (see README), so every example and test here works against
this synthetic catchment instead.

`build_raw_scenario` returns plain NumPy/Python values with no
`rri_torch` dependency; `build_synthetic_basin` (below) and
`reference.build_synthetic_basin_ref` (in `tests/test_reference_agreement.py`)
both build on top of it, so the differentiable model and the independent
loop-based reference are guaranteed to run the exact same catchment and
parameters.
"""

from __future__ import annotations

import numpy as np
import torch

DTYPE = torch.float64


def build_raw_scenario(
    ny: int = 12,
    nx: int = 7,
    dx: float = 50.0,
    dy: float = 50.0,
    along_slope: float = 0.01,
    cross_slope: float = 0.02,
    base_elev: float = 20.0,
    soildepth: float = 0.5,
    gammaa: float = 0.4,
    dm_frac: float = 0.3,
    ns_slope: float = 0.3,
    ka: float = 1.0e-5,
    beta: float = 5.0,
    min_wc4latflow: float = 0.0,
    river_width: float = 3.0,
    river_depth: float = 1.0,
    river_height: float = 0.0,
    ns_river: float = 0.03,
    ksv: float = 1.0e-6,
    faif: float = 0.1,
    infilt_limit: float = -1.0,
    outlet_slope: float = 1e-3,
) -> dict:
    center = nx // 2
    i_idx, j_idx = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")

    zb = base_elev - along_slope * i_idx * dy + cross_slope * np.abs(j_idx - center) * dx

    domain = np.ones((ny, nx), dtype=bool)
    riv_mask = np.zeros((ny, nx), dtype=bool)
    riv_mask[:, center] = True

    # D8 direction: every river cell flows south (code 4); irrelevant for
    # non-river cells since slope routing here is the full diffusive-wave
    # scheme, not the direction-based kinematic mode.
    dir_grid = np.full((ny, nx), 4, dtype=np.int64)
    len_riv = np.full((ny, nx), dy, dtype=np.float64)

    da = soildepth * gammaa
    dm = dm_frac * da

    return dict(
        ny=ny, nx=nx, dx=dx, dy=dy, center=center,
        zb=zb, domain=domain, riv_mask=riv_mask, dir_grid=dir_grid, len_riv=len_riv,
        soildepth=soildepth, gammaa=gammaa, da=da, dm=dm,
        ns_slope=ns_slope, ka=ka, beta=beta, min_wc4latflow=min_wc4latflow,
        river_width=river_width, river_depth=river_depth, river_height=river_height,
        ns_river=ns_river, ksv=ksv, faif=faif, infilt_limit=infilt_limit,
        outlet_slope=outlet_slope,
    )


def build_synthetic_basin(**kwargs):
    """Returns (grid, params, riv_outlet_index) for `rri_torch`."""
    from rri_torch.geometry import Grid
    from rri_torch.model import RRIParams
    from rri_torch.river import RiverParams
    from rri_torch.slope import SlopeParams

    raw = build_raw_scenario(**kwargs)
    ny, nx = raw["ny"], raw["nx"]

    grid = Grid.build(
        domain=raw["domain"], zb=raw["zb"], dir_grid=raw["dir_grid"],
        riv_mask=raw["riv_mask"], len_riv=raw["len_riv"],
        dx=raw["dx"], dy=raw["dy"], dtype=DTYPE,
    )
    n_riv = grid.n_riv
    outlet_index = int(np.argmax(grid.riv_i.numpy()))  # southernmost river cell

    def full(shape, value):
        return torch.full(shape, value, dtype=DTYPE)

    soildepth = full((ny, nx), raw["soildepth"])
    gammaa = full((ny, nx), raw["gammaa"])

    slope_params = SlopeParams(
        ns=full((ny, nx), raw["ns_slope"]),
        ka=full((ny, nx), raw["ka"]),
        beta=full((ny, nx), raw["beta"]),
        da=full((ny, nx), raw["da"]),
        dm=full((ny, nx), raw["dm"]),
        soildepth=soildepth,
        gammaa=gammaa,
        min_wc4latflow=full((ny, nx), raw["min_wc4latflow"]),
    )

    river_params = RiverParams(
        width=full((n_riv,), raw["river_width"]),
        depth=full((n_riv,), raw["river_depth"]),
        ns_river=torch.tensor(raw["ns_river"], dtype=DTYPE),
        height=full((n_riv,), raw["river_height"]),
        outlet_slope=raw["outlet_slope"],
    )

    params = RRIParams(
        slope=slope_params,
        river=river_params,
        ksv=full((ny, nx), raw["ksv"]),
        faif=full((ny, nx), raw["faif"]),
        infilt_limit=full((ny, nx), raw["infilt_limit"]),
        evp_switch=0,
    )

    return grid, params, outlet_index


def make_rain_pulse(
    T: int,
    ny: int,
    nx: int,
    peak_mm_per_hr: float = 30.0,
    start: int = 2,
    duration: int = 6,
) -> torch.Tensor:
    """A simple block-shaped rainfall pulse, uniform over the grid,
    returned as a rate in m/s for each of `T` steps."""
    rate = peak_mm_per_hr * 1e-3 / 3600.0
    rain = torch.zeros((T, ny, nx), dtype=DTYPE)
    rain[start:start + duration] = rate
    return rain
