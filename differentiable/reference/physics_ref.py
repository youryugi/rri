"""Loop-based (not vectorized) ports of the RHS/exchange/vertical-process
subroutines, mirroring the Fortran control flow cell-by-cell. Companion to
`hydraulics_ref.py`; together these are an independent cross-check of
`rri_torch`'s vectorized tensor implementation.
"""

from __future__ import annotations

import math

import numpy as np

from .hydraulics_ref import h2lev, hq_slope, hq_river
from .geometry_ref import RefGrid

_G = 9.81
_MU1 = (2.0 / 3.0) ** 1.5
_MU2 = 0.35
_MU3 = 0.91
N_CORRECTION_MAX = 10  # matches the Fortran `do count = 1, 10`


def slope_rhs_ref(hs: np.ndarray, grid: RefGrid, sp: dict, rain: np.ndarray) -> np.ndarray:
    """RRI_Slope.f90: `funcs` + `qs_calc_do`, restricted to eight_dir=1."""
    ny, nx = grid.ny, grid.nx
    lev = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            if grid.domain[i, j]:
                lev[i, j] = h2lev(hs[i, j], sp["soildepth"][i, j], sp["gammaa"][i, j])

    qedge = {}
    for i in range(ny):
        for j in range(nx):
            if not grid.domain[i, j]:
                continue
            for edge, ii, jj, distance, length in grid.neighbors_slope(i, j):
                dh = ((grid.zb[i, j] + lev[i, j]) - (grid.zb[ii, jj] + lev[ii, jj])) / distance
                if dh >= 0.0:
                    hw = hs[i, j]
                    if grid.zb[i, j] < grid.zb[ii, jj]:
                        hw = max(0.0, grid.zb[i, j] + hs[i, j] - grid.zb[ii, jj])
                    ka_p = sp["ka"][i, j]
                    min_hs = sp["soildepth"][i, j] * sp["gammaa"][i, j] * sp["min_wc4latflow"][i, j] if ka_p > 0.0 else 0.0
                    q = hq_slope(sp["ns"][i, j], ka_p, sp["da"][i, j], sp["dm"][i, j], sp["beta"][i, j],
                                 hw, dh, length, grid.area, min_hs)
                else:
                    hw = hs[ii, jj]
                    if grid.zb[ii, jj] < grid.zb[i, j]:
                        hw = max(0.0, grid.zb[ii, jj] + hs[ii, jj] - grid.zb[i, j])
                    ka_n = sp["ka"][ii, jj]
                    min_hs = sp["soildepth"][ii, jj] * sp["gammaa"][ii, jj] * sp["min_wc4latflow"][ii, jj] if ka_n > 0.0 else 0.0
                    q = -hq_slope(sp["ns"][ii, jj], ka_n, sp["da"][ii, jj], sp["dm"][ii, jj], sp["beta"][ii, jj],
                                  hw, dh, length, grid.area, min_hs)
                qedge[(i, j, edge)] = q

    fs = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            if not grid.domain[i, j]:
                continue
            outflow = sum(qedge[(i, j, edge)] for edge, _, _, _, _ in grid.neighbors_slope(i, j))
            inflow = sum(qedge[(ii, jj, edge)] for edge, ii, jj in grid.inflow_neighbors_slope(i, j))
            fs[i, j] = rain[i, j] - outflow + inflow
    return fs


def river_rhs_ref(hr: dict, grid: RefGrid, rp: dict) -> tuple[dict, dict]:
    """RRI_Riv.f90: `funcr` + `qr_calc`, rectangular channel only.

    `hr`, `rp["width"]`, `rp["depth"]`, `rp["ns_river"]` are dicts keyed by
    (i, j) river-cell coordinates. Returns (dhr_dt, qr) as matching dicts.
    """
    qr = {}
    for (i, j), down in grid.down_riv.items():
        zb_p = grid.zb[i, j] - rp["depth"][(i, j)]
        hr_p = hr[(i, j)]
        distance = grid.dis_riv[(i, j)]
        if down is not None:
            ii, jj = down
            zb_n = grid.zb[ii, jj] - rp["depth"][(ii, jj)]
            hr_n = hr[(ii, jj)]
        else:
            zb_n = zb_p - distance * rp["outlet_slope"]
            hr_n = 0.0

        dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance
        if dh >= 0.0:
            hw = hr_p
            if zb_p < zb_n:
                hw = max(0.0, zb_p + hr_p - zb_n)
            q = hq_river(hw, dh, rp["width"][(i, j)], rp["ns_river"])
        else:
            hw = hr_n
            if zb_n < zb_p:
                hw = max(0.0, zb_n + hr_n - zb_p)
            q = -hq_river(hw, dh, rp["width"][(i, j)], rp["ns_river"])
        qr[(i, j)] = q

    inflow = {k: 0.0 for k in grid.down_riv}
    for (i, j), down in grid.down_riv.items():
        if down is not None:
            inflow[down] += qr[(i, j)]

    dhr_dt = {}
    for (i, j) in grid.down_riv:
        channel_area = rp["width"][(i, j)] * grid.len_riv_grid[i, j]
        dhr_dt[(i, j)] = (-qr[(i, j)] + inflow[(i, j)]) / channel_area
    return dhr_dt, qr


def river_slope_exchange_ref(
    hr: dict, hs: np.ndarray, grid: RefGrid, rp: dict, dt: float,
) -> tuple[dict, np.ndarray, dict]:
    """RRI_RivSlo.f90: `funcrs_dt`, including the real `do count = 1, 10
    ... exit` overshoot-correction loop (this reference doesn't need to
    keep the computational graph static, so there's no reason to unroll
    it to a fixed count the way `rri_torch.exchange` does)."""
    hs = hs.copy()
    hr = dict(hr)
    qrs = {k: 0.0 for k in grid.down_riv}

    for (i, j) in grid.down_riv:
        height = max(rp["height"][(i, j)], 0.0)
        width = rp["width"][(i, j)]
        len_riv = grid.len_riv_grid[i, j]
        depth = rp["depth"][(i, j)]

        hs_top = hs[i, j]
        hr_top = hr[(i, j)] - depth

        if (height == 0.0 and hr_top < 0.0) or (height > 0.0 and hr_top < 0.0 and hs_top <= height):
            # case (a): slope -> river, free fall
            hrs = _MU1 * hs_top * math.sqrt(_G * max(hs_top, 0.0)) * dt * len_riv / grid.area
            hrs = min(hrs, hs_top)
            hs[i, j] -= hrs
            hr[(i, j)] += hrs * grid.area / (width * len_riv)
            qrs[(i, j)] = hrs

            hs_top = hs[i, j]
            hr_top = hr[(i, j)] - depth
            if hr_top >= -1e-5 and hr_top > hs_top:
                for _ in range(N_CORRECTION_MAX):
                    ar = len_riv * width / grid.area
                    hrs2 = (hs_top - hr_top) / (1.0 + 1.0 / ar)
                    hs[i, j] -= hrs2
                    hr[(i, j)] += hrs2 * grid.area / (width * len_riv)
                    qrs[(i, j)] += hrs2
                    if abs(hs[i, j] - (hr[(i, j)] - depth)) < 1e-5:
                        break
                    hs_top = hs[i, j]
                    hr_top = hr[(i, j)] - depth
                hr[(i, j)] = hs[i, j] + depth

        elif height > 0.0 and hs_top <= height and hr_top <= height and hr_top >= 0.0:
            # case (b): no exchange
            qrs[(i, j)] = 0.0

        elif hs_top <= hr_top and hr_top >= height:
            # case (c): river -> slope, overtopping
            h1 = hr_top - height
            h2 = hs_top - height
            if h1 <= 0.0:
                hrs = 0.0
            elif h2 / h1 <= 2.0 / 3.0:
                hrs = -_MU2 * h1 * math.sqrt(2.0 * _G * max(h1, 0.0)) * dt * len_riv / grid.area
            else:
                hrs = -_MU3 * h2 * math.sqrt(2.0 * _G * max(h1 - h2, 0.0)) * dt * len_riv / grid.area
            ar = len_riv * width / grid.area
            if abs(hrs) > abs(-(hr_top - height) * ar):
                hrs = -(hr_top - height) * ar
            qrs[(i, j)] = hrs

            hs[i, j] -= hrs
            hr[(i, j)] += hrs * grid.area / (width * len_riv)

            hs_top = hs[i, j]
            hr_top = hr[(i, j)] - depth
            if hs_top > hr_top:
                for _ in range(N_CORRECTION_MAX):
                    ar = len_riv * width / grid.area
                    hrs2 = (hs_top - hr_top) / (1.0 + 1.0 / ar)
                    hs[i, j] -= hrs2
                    hr[(i, j)] += hrs2 * grid.area / (width * len_riv)
                    qrs[(i, j)] += hrs2
                    if abs(hs[i, j] - (hr[(i, j)] - depth)) < 1e-5:
                        break
                    hs_top = hs[i, j]
                    hr_top = hr[(i, j)] - depth
                hr[(i, j)] = hs[i, j] + depth

        elif hs_top >= hr_top and hs_top >= height:
            # case (d): slope -> river, overtopping
            h1 = hs_top - height
            h2 = hr_top - height
            if h1 <= 0.0:
                hrs = 0.0
            elif h2 / h1 <= 2.0 / 3.0:
                hrs = _MU2 * h1 * math.sqrt(2.0 * _G * max(h1, 0.0)) * dt * len_riv / grid.area
            else:
                hrs = _MU3 * h2 * math.sqrt(2.0 * _G * max(h1 - h2, 0.0)) * dt * len_riv / grid.area
            if hrs > (hs_top - height):
                hrs = hs[i, j] - height
            qrs[(i, j)] = hrs

            hs[i, j] -= hrs
            hr[(i, j)] += hrs * grid.area / (width * len_riv)

            hs_top = hs[i, j]
            hr_top = hr[(i, j)] - depth
            if hr_top >= -1e-5 and hr_top > hs_top:
                for _ in range(N_CORRECTION_MAX):
                    ar = len_riv * width / grid.area
                    hrs2 = (hs_top - hr_top) / (1.0 + 1.0 / ar)
                    hs[i, j] -= hrs2
                    hr[(i, j)] += hrs2 * grid.area / (width * len_riv)
                    qrs[(i, j)] += hrs2
                    if abs(hs[i, j] - (hr[(i, j)] - depth)) < 1e-5:
                        break
                    hs_top = hs[i, j]
                    hr_top = hr[(i, j)] - depth
                hr[(i, j)] = hs[i, j] + depth
        else:
            qrs[(i, j)] = 0.0

    return hr, hs, qrs


def infiltration_ref(hs, gampt_ff, ksv, faif, gammaa, infilt_limit, dt):
    gampt_ff_temp = gampt_ff if gampt_ff > 0.01 else 0.01
    f = ksv * (1.0 + faif * gammaa / gampt_ff_temp)
    if f >= hs / dt:
        f = hs / dt
    if infilt_limit >= 0.0 and gampt_ff >= infilt_limit:
        f = 0.0
    gampt_ff_new = gampt_ff + f * dt
    hs_new = hs - f * dt
    if hs_new <= 0.0:
        hs_new = 0.0
    return hs_new, gampt_ff_new, f


def evaporation_ref(hs, gampt_ff, pet, dt, evp_switch):
    if evp_switch == 1:
        aevp = min(pet, (hs + gampt_ff) / dt)
    else:
        aevp = min(pet, hs / dt)
    hs_new = hs - aevp * dt
    if hs_new < 0.0:
        if evp_switch == 1:
            gampt_ff = gampt_ff + hs_new
        hs_new = 0.0
        if gampt_ff < 0.0:
            gampt_ff = 0.0
    return hs_new, gampt_ff, aevp
