"""Map surveyed Ishikari cross sections onto the RRI river grid.

The survey-derived table is produced by
``real_data/ishikari/build_channel_geometry.py``.  On covered mainstem
cells we use the surveyed bankfull depth and the width of an equivalent
rectangle that preserves bankfull cross-sectional area.  Power laws fitted
to those same sections provide a geometry fallback for tributaries and the
mainstem upstream of the survey coverage.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

REALDATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "real_data", "ishikari",
)
OUTLET_ROW_LOCAL, OUTLET_COL_LOCAL = 1663 - 1497, 10140 - 10096
OUTLET_KP = 26.6  # river.go.jp: Ishikari-Ohashi is 26.60km from the mouth

CODE_TO_DELTA = {1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
                 16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1)}
OPPOSITE_OFFSETS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def trace_mainstem(
    riv_mask: np.ndarray,
    dirarr: np.ndarray,
    accarr: np.ndarray,
    dx: float,
    dy: float,
):
    ny, nx = riv_mask.shape

    def upstream_neighbors(r, c):
        out = []
        for dr, dc in OPPOSITE_OFFSETS:
            nr, nc = r + dr, c + dc
            if nr < 0 or nr >= ny or nc < 0 or nc >= nx or not riv_mask[nr, nc]:
                continue
            delta = CODE_TO_DELTA.get(int(dirarr[nr, nc]))
            if delta is not None and (nr + delta[0], nc + delta[1]) == (r, c):
                out.append((nr, nc))
        return out

    r, c = OUTLET_ROW_LOCAL, OUTLET_COL_LOCAL
    cells, cum_km = [(r, c)], [0.0]
    while True:
        ups = upstream_neighbors(r, c)
        if not ups:
            break
        ups.sort(key=lambda rc: accarr[rc[0], rc[1]], reverse=True)
        nr, nc = ups[0]
        diagonal = (dx ** 2 + dy ** 2) ** 0.5
        step_km = (diagonal if (nr != r and nc != c) else (dy if nr != r else dx)) / 1000.0
        cum_km.append(cum_km[-1] + step_km)
        cells.append((nr, nc))
        r, c = nr, nc
    return cells, cum_km


def build_real_geometry_fields(
    ny: int,
    nx: int,
    riv_mask: np.ndarray,
    dirarr: np.ndarray,
    accarr: np.ndarray,
    dx: float,
    dy: float,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return per-grid-cell ``(width, depth, info)`` arrays."""
    survey = pd.read_csv(f"{REALDATA_DIR}/kp_channel_geometry.csv").sort_values("kp")
    kp_arr = survey["kp"].to_numpy()
    surveyed_width = survey["equivalent_width_m"].to_numpy()
    surveyed_depth = survey["depth_m"].to_numpy()
    mainstem, cum_km = trace_mainstem(riv_mask, dirarr, accarr, dx, dy)
    absolute_kp = OUTLET_KP + np.array(cum_km)
    kp_min, kp_max = kp_arr.min(), kp_arr.max()
    covered = (absolute_kp >= kp_min) & (absolute_kp <= kp_max)
    if covered.sum() < 2:
        raise ValueError("Too few surveyed mainstem cells to fit fallback geometry")

    mainstem_area_km2 = np.array([accarr[r, c] for r, c in mainstem]) * dx * dy * 1e-6
    width_on_grid = np.interp(absolute_kp[covered], kp_arr, surveyed_width)
    depth_on_grid = np.interp(absolute_kp[covered], kp_arr, surveyed_depth)

    def fit_power_law(values: np.ndarray) -> tuple[float, float, float]:
        log_area = np.log(mainstem_area_km2[covered])
        log_values = np.log(values)
        exponent, log_coefficient = np.polyfit(log_area, log_values, 1)
        predicted = log_coefficient + exponent * log_area
        residual = np.sum((log_values - predicted) ** 2)
        total = np.sum((log_values - log_values.mean()) ** 2)
        r2 = 1.0 - residual / total
        return float(np.exp(log_coefficient)), float(exponent), float(r2)

    width_c, width_s, width_r2 = fit_power_law(width_on_grid)
    depth_c, depth_s, depth_r2 = fit_power_law(depth_on_grid)
    drainage_km2 = np.maximum(accarr * dx * dy * 1e-6, 1e-12)
    width_field = np.zeros((ny, nx), dtype=float)
    depth_field = np.zeros((ny, nx), dtype=float)
    width_field[riv_mask] = width_c * drainage_km2[riv_mask] ** width_s
    depth_field[riv_mask] = depth_c * drainage_km2[riv_mask] ** depth_s

    for index, ((r, c), is_covered) in enumerate(zip(mainstem, covered)):
        if is_covered:
            width_field[r, c] = float(np.interp(absolute_kp[index], kp_arr, surveyed_width))
            depth_field[r, c] = float(np.interp(absolute_kp[index], kp_arr, surveyed_depth))

    info = {
        "mainstem_cells": len(mainstem),
        "mainstem_length_km": cum_km[-1],
        "real_data_cells": int(covered.sum()),
        "formula_fallback_cells": int((~covered).sum()),
        "kp_coverage": (kp_min, kp_max),
        "width_power_law": (width_c, width_s, width_r2),
        "depth_power_law": (depth_c, depth_s, depth_r2),
    }
    return width_field, depth_field, info
