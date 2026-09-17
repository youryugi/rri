"""Loader for the real Solo River basin project (`RRI-CUI/Project/solo30s`),
the basin used in RRI's own official tutorial (Sayama & Iwami, 2016).

Every parsing/geometry rule here was cross-checked against the actual
Fortran source (`RRI_1_4_2_7_GUI_Beta/RRI-CUI/source/1.4.2.7/*.f90`), not
guessed from the file formats alone -- see the docstring of each function
for the exact source location it mirrors. This file has no `rri_torch`
import in its parsing helpers; only `build_solo_scenario` at the bottom
wires the parsed arrays into `rri_torch`'s `Grid`/`RRIParams`.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# ESRI ASCII grid (adem.txt / acc_mod.txt / dir_mod.txt)
# ---------------------------------------------------------------------------

@dataclass
class EsriGrid:
    ncols: int
    nrows: int
    xllcorner: float
    yllcorner: float
    cellsize: float
    nodata: float
    data: np.ndarray  # (nrows, ncols), row 0 = north (matches RRI's (i, j) with i=1 at north)


def read_esri_ascii(path: str) -> EsriGrid:
    """Parses the 6-line ESRI ASCII header + comma-separated rows used by
    every file under solo30s/topo/. Mirrors the plain `read(10,*) ctemp, x`
    header reads and the `read_gis_int`/`read_gis_real` subroutines in
    RRI_Sub.f90 (row 1 of the file = northernmost row = i=1)."""
    with open(path) as f:
        header = {}
        for _ in range(6):
            key, val = f.readline().split(None, 1)
            header[key.strip().lower()] = val.strip()
        ncols = int(header["ncols"])
        nrows = int(header["nrows"])
        xllcorner = float(header["xllcorner"])
        yllcorner = float(header["yllcorner"])
        cellsize = float(header["cellsize"])
        nodata = float(header["nodata_value"])

        data = np.empty((nrows, ncols), dtype=np.float64)
        for i in range(nrows):
            line = f.readline()
            vals = line.replace(",", " ").split()
            data[i, :] = [float(v) for v in vals]

    return EsriGrid(ncols=ncols, nrows=nrows, xllcorner=xllcorner,
                     yllcorner=yllcorner, cellsize=cellsize, nodata=nodata,
                     data=data)


# ---------------------------------------------------------------------------
# Hubeny geodesic distance (lat-lon dx/dy), mirrors `hubeny_sub` in
# RRI_Sub.f90 exactly (GRS80 ellipsoid constants included).
# ---------------------------------------------------------------------------

_A = 6378137.000        # GRS80 semi-major axis [m]
_B = 6356752.314        # GRS80 semi-minor axis [m]


def hubeny_distance(x1_deg: float, y1_deg: float, x2_deg: float, y2_deg: float) -> float:
    """Geodesic distance in metres between two lon/lat points, exactly
    reproducing `hubeny_sub` (RRI_Sub.f90, line ~499)."""
    x1, y1 = math.radians(x1_deg), math.radians(y1_deg)
    x2, y2 = math.radians(x2_deg), math.radians(y2_deg)
    dy = y1 - y2
    dx = x1 - x2
    mu = (y1 + y2) / 2.0
    e = math.sqrt((_A ** 2 - _B ** 2) / _A ** 2)
    w = math.sqrt(1.0 - e ** 2 * math.sin(mu) ** 2)
    n = _A / w
    m = _A * (1.0 - e ** 2) / w ** 3
    return math.sqrt((dy * m) ** 2 + (dx * n * math.cos(mu)) ** 2)


def compute_dxdy(nx: int, ny: int, xllcorner: float, yllcorner: float,
                  cellsize: float, utm: bool) -> tuple[float, float]:
    """Single scalar (dx, dy) in metres for the whole domain, mirroring
    RRI.f90's "STEP 2: CALC PREPARATION" (lines ~117-152): for utm grids,
    dx=dy=cellsize; for lat-lon grids, average the Hubeny distance of the
    south/north edges for dx and the west/east edges for dy."""
    if utm:
        return cellsize, cellsize
    x0, y0 = xllcorner, yllcorner
    x1, y1 = xllcorner + nx * cellsize, yllcorner + ny * cellsize
    d1 = hubeny_distance(x0, y0, x1, y0)          # south edge
    d2 = hubeny_distance(x0, y1, x1, y1)          # north edge
    d3 = hubeny_distance(x0, y0, x0, y1)          # west edge
    d4 = hubeny_distance(x1, y0, x1, y1)          # east edge
    dx = (d1 + d2) / 2.0 / nx
    dy = (d3 + d4) / 2.0 / ny
    return dx, dy


# ---------------------------------------------------------------------------
# rain.dat: repeating blocks of `t_seconds nx_rain ny_rain` + ny_rain rows.
# Mirrors the read loop in RRI.f90 (lines ~500-535) including the
# grid-alignment index mapping and the step-function time lookup.
# ---------------------------------------------------------------------------

@dataclass
class RainData:
    t_block: np.ndarray     # (n_blocks,) block end-times [s], t_block[0] == 0
    values: np.ndarray      # (n_blocks, ny_rain, nx_rain) [m/s]
    nx_rain: int
    ny_rain: int
    xllcorner: float
    yllcorner: float
    cellsize_x: float
    cellsize_y: float

    def rain_index_map(self, nx: int, ny: int, xllcorner: float,
                        yllcorner: float, cellsize: float) -> tuple[np.ndarray, np.ndarray]:
        """1-indexed-in-Fortran, returned here as 0-indexed numpy arrays:
        rain_j[j] (0-indexed j) and rain_i[i] (0-indexed i), the rain-grid
        cell each topo-grid cell reads from. Exactly mirrors lines 529-534
        of RRI.f90. Out-of-range entries are left as -1 (RRI skips them,
        leaving rain at 0 there)."""
        j = np.arange(nx)
        rain_j = np.floor((xllcorner + (j + 0.5) * cellsize - self.xllcorner) / self.cellsize_x).astype(np.int64)
        i = np.arange(ny)
        rain_i = self.ny_rain - 1 - np.floor(
            (yllcorner + (ny - 1 - i + 0.5) * cellsize - self.yllcorner) / self.cellsize_y
        ).astype(np.int64)
        rain_j[(rain_j < 0) | (rain_j >= self.nx_rain)] = -1
        rain_i[(rain_i < 0) | (rain_i >= self.ny_rain)] = -1
        return rain_i, rain_j

    def block_index_at(self, time_plus_ddt: float) -> int:
        """`itemp` from RRI.f90 line 752: the block whose *end* time is the
        first >= time_plus_ddt (blocks are right-continuous: block k's
        rainfall rate applies over `(t_block[k-1], t_block[k]]`)."""
        idx = int(np.searchsorted(self.t_block, time_plus_ddt, side="left"))
        return min(idx, len(self.t_block) - 1)


def read_rain_dat(path: str) -> RainData:
    blocks_t = []
    blocks_v = []
    nx_rain = ny_rain = None
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            t_str, nxr_str, nyr_str = line.split()
            t = float(t_str)
            nx_rain = int(nxr_str)
            ny_rain = int(nyr_str)
            grid = np.empty((ny_rain, nx_rain), dtype=np.float64)
            for i in range(ny_rain):
                grid[i, :] = [float(v) for v in f.readline().split()]
            blocks_t.append(t)
            blocks_v.append(grid)

    t_block = np.array(blocks_t, dtype=np.float64)
    values = np.stack(blocks_v, axis=0) / 3600.0 / 1000.0  # mm/h -> m/s
    return RainData(t_block=t_block, values=values, nx_rain=nx_rain, ny_rain=ny_rain,
                     xllcorner=None, yllcorner=None, cellsize_x=None, cellsize_y=None)


# ---------------------------------------------------------------------------
# obs/*.data and obs/location_solo_30s_all.txt
# ---------------------------------------------------------------------------

def read_obs_series(path: str, step_hours: Optional[float] = None) -> tuple[np.ndarray, np.ndarray]:
    """Plain `<time> <value>` pairs, one per line, e.g. `obs/disc_*.data`.

    `step_hours`: the first column is **not** simulation-hours -- it's
    the `calcHydro.f90` post-processing tool's output-step index
    (continuous-valued here, since observations don't land exactly on
    model output steps), where step index relates to hours as
    `hour = step_index * (lasth / outnum)` (`RRI-CUI/etc/calcHydro/calcHydro.f90`
    writes `t` -- a raw step counter, never converted to hours -- as the
    first column of its own `hydro_<name>.txt` output; `disc_*.data`'s
    first column is the same unit). For Solo (`lasth=360`, `outnum=96`)
    this is 3.75, and e.g. `disc_cepu.data`'s `0, 6.4, 12.8, ..., 96`
    becomes exactly `0, 24, 48, ..., 360` -- one observation per day of
    the 15-day event (see HANDOFF.md section 8.3 for the full derivation
    and NSE-based verification: reading the column as hours directly
    gives NSE=-2.32 against the official simulated hydrograph; converting
    with `step_hours=3.75` gives NSE=0.295). Pass the project's own
    `lasth/outnum` here -- default `None` returns the raw column
    unconverted (e.g. if the caller wants to match `calcHydro.f90`'s own
    output-step convention directly instead of real hours).
    """
    idx, val = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            a, b = line.split()
            idx.append(float(a))
            val.append(float(b))
    idx = np.array(idx)
    if step_hours is not None:
        idx = idx * step_hours
    return idx, np.array(val)


def read_locations(path: str) -> dict[str, tuple[int, int]]:
    """`name row col` per line, 1-indexed (i, j) as in RRI's convention."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, row, col = line.split()
            out[name] = (int(row), int(col))
    return out


def read_point_hydrograph(path: str) -> tuple[np.ndarray, np.ndarray]:
    """`hydro.txt` / `hydro_hr.txt`: pre-computed point time series at the
    single gauge in `location.txt`, written directly by RRI.f90 itself
    (search `hydro_file` in RRI.f90) -- NOT derived from the `out/qr_*.out`
    spatial dumps, and at hourly resolution regardless of `outnum`."""
    t, v = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            a, b = line.split()
            t.append(float(a))
            v.append(float(b))
    return np.array(t), np.array(v)


# ---------------------------------------------------------------------------
# Full scenario assembly
# ---------------------------------------------------------------------------

def build_solo_scenario(project_dir: str, dtype=None):
    """Builds (grid, params, outlet_riv_indices, rain_data) for the Solo
    River basin project, ready to hand to `rri_torch.model.RRIModel`.

    `outlet_riv_indices` is a dict {gauge_name: river-cell index into
    `grid`} for every station in `obs/location_solo_30s_all.txt`.

    Parameter choices here are read directly from `RRI_Input.txt` (see
    HANDOFF.md table); this project has ka=0 and ksv=0, so subsurface/
    infiltration terms are numerically inert -- only Manning surface flow
    (slope + river) is exercised, which is exactly what `rri_torch`
    implements.
    """
    import torch
    from rri_torch.geometry import Grid
    from rri_torch.model import RRIParams
    from rri_torch.river import RiverParams
    from rri_torch.slope import SlopeParams

    if dtype is None:
        dtype = torch.float64

    dem = read_esri_ascii(os.path.join(project_dir, "topo", "adem.txt"))
    acc = read_esri_ascii(os.path.join(project_dir, "topo", "acc_mod.txt"))
    dirg = read_esri_ascii(os.path.join(project_dir, "topo", "dir_mod.txt"))

    ny, nx = dem.nrows, dem.ncols
    zs = dem.data
    acc_arr = acc.data
    dir_arr = dirg.data.astype(np.int64)

    # domain / sink, mirrors RRI.f90 lines 234-244 exactly (zs > -100, not
    # `== nodata`, in case of near-nodata interpolation artifacts).
    domain = zs > -100.0
    sink = domain & ((dir_arr == 0) | (dir_arr == -1))

    dx, dy = compute_dxdy(nx, ny, dem.xllcorner, dem.yllcorner, dem.cellsize, utm=False)
    area = dx * dy
    length = math.sqrt(dx * dy)

    # RRI_Input.txt parameters for this project (see HANDOFF.md table).
    riv_thresh = 100.0
    width_param_c, width_param_s = 5.000, 0.350
    depth_param_c, depth_param_s = 0.950, 0.200
    height_param, height_limit_param = 0.000, 20.0
    ns_slope = 0.400
    ns_river = 0.030
    soildepth_val = 1.000
    gammaa_val = 0.475
    ka_val = 0.000
    ksv_val = 0.000
    faif_val = 0.316
    beta_val = 8.000

    riv_mask = domain & (acc_arr > riv_thresh)
    width = np.zeros((ny, nx))
    depth = np.zeros((ny, nx))
    height = np.zeros((ny, nx))
    drainage_km2 = acc_arr * dx * dy * 1e-6
    width[riv_mask] = width_param_c * drainage_km2[riv_mask] ** width_param_s
    depth[riv_mask] = depth_param_c * drainage_km2[riv_mask] ** depth_param_s
    height[riv_mask & (acc_arr > height_limit_param)] = height_param
    len_riv_grid = np.where(riv_mask, length, 0.0)

    grid = Grid.build(
        domain=domain, zb=zs, dir_grid=dir_arr,
        riv_mask=riv_mask, sink=sink, len_riv=len_riv_grid,
        dx=dx, dy=dy, dtype=dtype,
    )
    n_riv = grid.n_riv

    def full(shape, value):
        return torch.full(shape, value, dtype=dtype)

    soildepth = full((ny, nx), soildepth_val)
    gammaa = full((ny, nx), gammaa_val)
    # da = soildepth*gammaa only if ka>0 else 0 (RRI_Read.f90 line 300);
    # here ka=0 so da=dm=0 -- lateral unsaturated flow term is inert.
    da = torch.zeros((ny, nx), dtype=dtype)
    dm = torch.zeros((ny, nx), dtype=dtype)

    slope_params = SlopeParams(
        ns=full((ny, nx), ns_slope),
        ka=full((ny, nx), ka_val),
        beta=full((ny, nx), beta_val),
        da=da, dm=dm,
        soildepth=soildepth, gammaa=gammaa,
        min_wc4latflow=full((ny, nx), 0.0),
    )

    river_params = RiverParams(
        width=grid.gather_to_riv(torch.as_tensor(width, dtype=dtype)),
        depth=grid.gather_to_riv(torch.as_tensor(depth, dtype=dtype)),
        ns_river=torch.tensor(ns_river, dtype=dtype),
        height=grid.gather_to_riv(torch.as_tensor(height, dtype=dtype)),
        outlet_slope=1e-3,
    )

    params = RRIParams(
        slope=slope_params,
        river=river_params,
        ksv=full((ny, nx), ksv_val),
        faif=full((ny, nx), faif_val),
        infilt_limit=full((ny, nx), -1.0),
        evp_switch=0,
    )

    locations = read_locations(os.path.join(project_dir, "obs", "location_solo_30s_all.txt"))
    outlet_riv_indices = {}
    riv_index_map = grid.riv_index_map.numpy()
    for name, (row, col) in locations.items():
        i, j = row - 1, col - 1  # RRI is 1-indexed
        idx = riv_index_map[i, j]
        if idx >= 0:
            outlet_riv_indices[name] = int(idx)

    rain = read_rain_dat(os.path.join(project_dir, "rain", "rain.dat"))
    rain.xllcorner = dem.xllcorner
    rain.yllcorner = dem.yllcorner
    rain.cellsize_x = dem.cellsize
    rain.cellsize_y = dem.cellsize

    return grid, params, outlet_riv_indices, rain


def build_rain_sequence(rain: RainData, nx: int, ny: int, xllcorner: float,
                         yllcorner: float, cellsize: float, dt: float, n_steps: int,
                         dtype=None) -> "torch.Tensor":
    """Materializes a (n_steps, ny, nx) rain-rate tensor [m/s], sampling
    RRI's step-function rain (see `RainData.block_index_at`) at the end
    time of each step, matching RRI's own `time + ddt` bin lookup when
    `dt` divides evenly into the block spacing (true for Solo: 600s steps,
    86400s blocks)."""
    import torch
    if dtype is None:
        dtype = torch.float64

    rain_i, rain_j = rain.rain_index_map(nx, ny, xllcorner, yllcorner, cellsize)
    valid_i = rain_i >= 0
    valid_j = rain_j >= 0

    out = np.zeros((n_steps, ny, nx), dtype=np.float64)
    ii = np.where(valid_i)[0]
    jj = np.where(valid_j)[0]
    if len(ii) and len(jj):
        ri = rain_i[ii]
        rj = rain_j[jj]
        for t in range(n_steps):
            block = rain.block_index_at((t + 1) * dt)
            sub = rain.values[block][np.ix_(ri, rj)]
            out[np.ix_([t], ii, jj)] = sub[None, :, :]

    return torch.as_tensor(out, dtype=dtype)
