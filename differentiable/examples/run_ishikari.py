"""Event simulation and numerical diagnostics for the Ishikari River.

No RRI_Input.txt / rain.dat files are written: the script builds the
Grid/RRIParams directly in Python, maps official surveyed cross sections
onto the mainstem, and uses fitted width/depth fallbacks elsewhere.  Daily
basin-averaged MERV-Jp precipitation is broadcast uniformly, consistent
with how that source forcing was derived.

There is no official-Fortran reference for this basin (no compiled
example ships for it, unlike Solo) -- this can only ever be validated
against the CSV's real "Obs flow" column, not a second "vs official"
axis.

Run with:  python examples/run_ishikari.py --device cuda
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from real_data.solo_river import read_esri_ascii, compute_dxdy
from rri_torch.geometry import Grid
from rri_torch.slope import SlopeParams
from rri_torch.river import RiverParams
from rri_torch.model import RRIModel, RRIParams

PROJECT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "..", "RRI_1_4_2_7_GUI_Beta", "RRI-CUI", "Project", "ishikari")
CATCHMENT_AREA_KM2 = 12697.0  # official GRDC/MLIT catchment area at Ishikari-Ohashi


def build_ishikari_scenario(dtype=torch.float64):
    dem = read_esri_ascii(os.path.join(PROJECT_DIR, "topo", "adem.txt"))
    acc = read_esri_ascii(os.path.join(PROJECT_DIR, "topo", "acc_mod.txt"))
    dirg = read_esri_ascii(os.path.join(PROJECT_DIR, "topo", "dir_mod.txt"))

    ny, nx = dem.nrows, dem.ncols
    zs = dem.data
    acc_arr = acc.data
    dir_arr = dirg.data.astype(np.int64)

    domain = zs > -100.0
    sink = domain & ((dir_arr == 0) | (dir_arr == -1))

    dx, dy = compute_dxdy(nx, ny, dem.xllcorner, dem.yllcorner, dem.cellsize, utm=False)
    length = math.sqrt(dx * dy)
    print(f"dx={dx:.1f}m dy={dy:.1f}m cell_area={dx*dy/1e6:.4f}km2 active_cells={int(domain.sum())}")

    # The channel-geometry helper fits width/depth power laws to official
    # Ishikari cross sections for fallback use, then applies the surveyed
    # equivalent rectangles directly on the covered mainstem cells.
    riv_thresh = 100.0
    ns_slope, ns_river = 0.400, 0.030

    riv_mask = domain & (acc_arr > riv_thresh)
    print(f"river cells: {int(riv_mask.sum())} (threshold acc>{riv_thresh})")
    from examples.ishikari_real_geometry import build_real_geometry_fields
    width, depth, geom_info = build_real_geometry_fields(
        ny, nx, riv_mask, dir_arr, acc_arr, dx=dx, dy=dy,
    )
    wc, ws, wr2 = geom_info["width_power_law"]
    dc, ds, dr2 = geom_info["depth_power_law"]
    print(f"survey geometry: {geom_info['mainstem_cells']} mainstem cells, "
          f"{geom_info['mainstem_length_km']:.1f}km traced; "
          f"{geom_info['real_data_cells']} cells got surveyed equivalent rectangles "
          f"(KP {geom_info['kp_coverage'][0]:.1f}-{geom_info['kp_coverage'][1]:.1f}), "
          f"{geom_info['formula_fallback_cells']} mainstem cells use fitted fallback")
    print(f"fitted fallback: width={wc:.3f}*A^{ws:.3f} (log-R2={wr2:.3f}), "
          f"depth={dc:.3f}*A^{ds:.3f} (log-R2={dr2:.3f}); A in km2")

    len_riv_grid = np.where(riv_mask, length, 0.0)

    grid = Grid.build(domain=domain, zb=zs, dir_grid=dir_arr, riv_mask=riv_mask,
                       sink=sink, len_riv=len_riv_grid, dx=dx, dy=dy, dtype=dtype)

    def full(shape, v):
        return torch.full(shape, v, dtype=dtype)

    slope_params = SlopeParams(
        ns=full((ny, nx), ns_slope), ka=full((ny, nx), 0.0), beta=full((ny, nx), 8.0),
        da=torch.zeros((ny, nx), dtype=dtype), dm=torch.zeros((ny, nx), dtype=dtype),
        soildepth=full((ny, nx), 1.0), gammaa=full((ny, nx), 0.475),
        min_wc4latflow=full((ny, nx), 0.0),
    )
    river_params = RiverParams(
        width=grid.gather_to_riv(torch.as_tensor(width, dtype=dtype)),
        depth=grid.gather_to_riv(torch.as_tensor(depth, dtype=dtype)),
        ns_river=torch.tensor(ns_river, dtype=dtype),
        height=grid.gather_to_riv(torch.zeros(ny, nx, dtype=dtype)),
        outlet_slope=1e-3,
    )
    params = RRIParams(slope=slope_params, river=river_params, ksv=full((ny, nx), 0.0),
                        faif=full((ny, nx), 0.316), infilt_limit=full((ny, nx), -1.0), evp_switch=0)

    # outlet: the snapped Ishikari-Ohashi pixel, in this clipped grid's
    # local row/col (see HANDOFF.md section 9 for the full-raster snap).
    r0_full, c0_full = 1497, 10096  # export_rri.py's clip window origin
    outlet_row_local = 1663 - r0_full
    outlet_col_local = 10140 - c0_full
    outlet_idx = int(grid.riv_index_map[outlet_row_local, outlet_col_local])
    assert outlet_idx >= 0, "outlet cell isn't marked as a river cell -- check riv_thresh"

    return grid, params, outlet_idx


def load_forcing(start_date: str, days: int):
    df = pd.read_csv(os.path.join(PROJECT_DIR, "obs", "varssim005_ver2_1.csv"))
    df["date"] = pd.to_datetime(df[["Year", "Month", "Day"]])
    df = df.set_index("date")
    window = df.loc[start_date:].iloc[:days]
    precip_mm_day = window["Precip"].to_numpy()
    obs_flow_mm_day = window["Obs flow"].to_numpy()
    dates = window.index
    obs_flow_m3s = obs_flow_mm_day * CATCHMENT_AREA_KM2 * 1000.0 / 86400.0
    return dates, precip_mm_day, obs_flow_m3s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2001-09-01")
    ap.add_argument("--days", type=int, default=20)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--adaptive", action="store_true", default=True)
    ap.add_argument("--eps", type=float, default=0.01)
    ap.add_argument("--ddt-min-slope", type=float, default=1.0)
    ap.add_argument("--ddt-min-river", type=float, default=0.1)
    args = ap.parse_args()

    device = torch.device(args.device)
    dtype = torch.float64

    grid, params, outlet_idx = build_ishikari_scenario(dtype=dtype)
    dates, precip_mm_day, obs_flow_m3s = load_forcing(args.start, args.days)
    print(f"forcing window: {dates[0].date()} to {dates[-1].date()}, "
          f"precip range {precip_mm_day.min():.1f}-{precip_mm_day.max():.1f} mm/day")

    grid = grid.to(device)
    import dataclasses
    def to_dev(o):
        if torch.is_tensor(o):
            return o.to(device)
        if dataclasses.is_dataclass(o):
            return dataclasses.replace(o, **{f.name: to_dev(getattr(o, f.name)) for f in dataclasses.fields(o)})
        return o
    params.slope = to_dev(params.slope)
    params.river = to_dev(params.river)
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)

    dt = 600.0
    steps_per_day = int(round(86400.0 / dt))
    T = args.days * steps_per_day
    rain_rate_m_s = precip_mm_day * 1e-3 / 86400.0  # mm/day -> m/s

    # Stream the forcing one outer step at a time.  Besides avoiding a
    # days*144*ny*nx rainfall allocation, this records the state maxima and
    # adaptive-step workload needed to distinguish a healthy flood response
    # from the pathological stiffness seen with the original unconditioned
    # DEM / borrowed channel geometry (HANDOFF.md section 10).
    model = RRIModel(
        grid, adaptive=True, eps=args.eps,
        ddt_min_slope=args.ddt_min_slope, ddt_min_river=args.ddt_min_river,
        track_qr_avg=True,
    )
    hs, gampt_ff, hr = model.initial_state(dtype=dtype)
    rain_t = torch.zeros((grid.ny, grid.nx), dtype=dtype, device=device)
    qr, hs_max, hr_max = [], [], []
    river_substeps, slope_substeps = [], []

    t0 = time.time()
    with torch.no_grad():
        for day in range(args.days):
            rain_t.zero_()
            rain_t[grid.domain] = rain_rate_m_s[day]
            day_start = time.time()
            for _ in range(steps_per_day):
                hs, gampt_ff, hr, diag = model.step(
                    hs, gampt_ff, hr, params, rain_t, None, dt,
                )
                q = diag.qr_avg if diag.qr_avg is not None else diag.qr
                qr.append(float(q[outlet_idx]))
                hs_max.append(float(hs[grid.domain].max()))
                hr_max.append(float(hr.max()))
                river_substeps.append(diag.river_stats.n_accepted)
                slope_substeps.append(diag.slope_stats.n_accepted)
            if device.type == "cuda":
                torch.cuda.synchronize()
            sl = slice(day * steps_per_day, (day + 1) * steps_per_day)
            print(
                f"day {day + 1:02d}/{args.days} {dates[day].date()} "
                f"P={precip_mm_day[day]:6.2f}mm "
                f"Qmean={np.mean(qr[sl]):8.1f} Qmax={np.max(qr[sl]):8.1f}m3/s "
                f"hs_max={np.max(hs_max[sl]):6.3f}m hr_max={np.max(hr_max[sl]):6.3f}m "
                f"substeps(riv/slo)={np.max(river_substeps[sl])}/{np.max(slope_substeps[sl])} "
                f"wall={time.time() - day_start:.1f}s",
                flush=True,
            )

    elapsed = time.time() - t0
    qr = np.asarray(qr)
    hs_max = np.asarray(hs_max)
    hr_max = np.asarray(hr_max)
    river_substeps = np.asarray(river_substeps)
    slope_substeps = np.asarray(slope_substeps)
    daily_qr = qr.reshape(args.days, steps_per_day).mean(axis=1)
    finite = np.isfinite(qr).all() and np.isfinite(hs_max).all() and np.isfinite(hr_max).all()
    print(f"forward sim done in {elapsed:.1f}s (eps={args.eps})")
    print(f"finite: {finite}, peak sim qr: {qr.max():.1f} m3/s at step {qr.argmax()}/{T}")
    print(f"state peaks: hs={hs_max.max():.3f}m, hr={hr_max.max():.3f}m; "
          f"max accepted substeps river={river_substeps.max()}, slope={slope_substeps.max()}")

    residual = daily_qr - obs_flow_m3s
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    denominator = float(np.sum((obs_flow_m3s - obs_flow_m3s.mean()) ** 2))
    nse = 1.0 - float(np.sum(residual ** 2)) / denominator if denominator > 0 else float("nan")
    print(f"daily comparison (cold-start, uncalibrated): RMSE={rmse:.1f}m3/s NSE={nse:.3f}")

    sim_daily_time = (np.arange(1, T+1) * dt) / 86400.0  # days since start
    obs_daily_time = np.arange(args.days) + 0.5  # obs "Obs flow" is a daily value, centered

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(sim_daily_time, qr, "-", color="#1f77b4", alpha=0.7, label="rri_torch (10-min avg)")
    ax.plot(obs_daily_time, daily_qr, "s-", color="#1f77b4", label="rri_torch (daily avg)")
    ax.plot(obs_daily_time, obs_flow_m3s, "o-", color="#d62728", label="Real observed (MERV-Jp)")
    ax2 = ax.twinx()
    ax2.bar(np.arange(args.days) + 0.5, precip_mm_day, width=0.8, color="gray", alpha=0.3, label="Precip (mm/day)")
    ax2.set_ylabel("Precipitation [mm/day]")
    ax2.invert_yaxis()
    ax.set_xlabel(f"Days since {dates[0].date()}")
    ax.set_ylabel(r"Discharge at Ishikari-Ohashi [m$^3$/s]")
    ax.set_title(
        f"Ishikari River: surveyed-geometry event simulation "
        f"({dates[0].date()} to {dates[-1].date()})"
    )
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=9)
    fig.tight_layout()

    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(out_dir, exist_ok=True)
    stem = f"ishikari_{dates[0].date()}_{args.days}d"
    figure_path = os.path.join(out_dir, f"fig_{stem}.png")
    data_path = os.path.join(out_dir, f"{stem}.npz")
    fig.savefig(figure_path, dpi=150)
    np.savez_compressed(
        data_path,
        dates=np.asarray(dates.astype(str), dtype="U10"), precip_mm_day=precip_mm_day,
        obs_flow_m3s=obs_flow_m3s, qr_10min_m3s=qr, qr_daily_m3s=daily_qr,
        hs_max_m=hs_max, hr_max_m=hr_max,
        river_accepted_substeps=river_substeps,
        slope_accepted_substeps=slope_substeps,
        elapsed_seconds=elapsed,
    )
    print(f"saved {figure_path}")
    print(f"saved {data_path}")


if __name__ == "__main__":
    main()
