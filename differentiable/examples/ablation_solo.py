"""Clean, controlled 2x2 ablation of the two formula fixes documented in
HANDOFF.md section 7.2 (hydraulic radius in `hq_river`, the "both banks"
factor of 2 in river<->slope exchange), run with IDENTICAL settings
(adaptive Cash-Karp integration, stage-accurate `qr_avg` tracking) for a
fair paper-table comparison -- unlike the historical runs during
debugging, which used inconsistent substep/qr_avg settings across fixes.

One `outlet_riv_index` list carries both the 6 named gauges AND the full
148-cell Cepu-to-outlet chain, so a single forward pass per variant
yields everything needed for both the hydrograph and the spatial-bias
figures -- no extra simulation runs.

Produces (under differentiable/results/):
  - ablation_table.csv / .md   : NSE/RMSE/bias/peak/runtime per variant
  - fig_hydrograph.png         : Cepu hydrograph, all 4 variants vs hydro.txt
  - fig_spatial_bias.png       : %diff along the chain, all 4 variants, at t=180h

Usage:
  python examples/ablation_solo.py --steps 2160 --device cuda
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import rri_torch.river as river_mod
import rri_torch.model as model_mod
from rri_torch.model import RRIModel
from real_data.solo_river import build_solo_scenario, build_rain_sequence, read_esri_ascii, read_point_hydrograph
from _legacy_formulas import legacy_hq_river, legacy_river_slope_exchange

PROJECT_DIR = "/home/yang/github/rri/RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s"
FORTRAN_OUT_DIR = "/tmp/claude-1004/-home-yang-github-rri/56b4d09e-3c12-4eb2-94cb-037713c8deea/scratchpad/solo30s_run/out"
DT = 600.0
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")

VARIANTS = {
    "both_legacy":  dict(river_fixed=False, exchange_fixed=False, label="Pre-fix (both legacy)"),
    "river_fixed":  dict(river_fixed=True,  exchange_fixed=False, label="Hydraulic-radius fix only"),
    "exchange_fixed": dict(river_fixed=False, exchange_fixed=True, label="Exchange x2 fix only"),
    "both_fixed":   dict(river_fixed=True,  exchange_fixed=True, label="Both fixes (current model)"),
}
COLORS = {"both_legacy": "#d62728", "river_fixed": "#ff7f0e", "exchange_fixed": "#9467bd", "both_fixed": "#1f77b4"}


def nash_sutcliffe(obs, sim):
    return 1.0 - np.sum((obs - sim) ** 2) / np.sum((obs - obs.mean()) ** 2)


def rmse(obs, sim):
    return float(np.sqrt(np.mean((obs - sim) ** 2)))


def fortran_nint(x):
    return math.floor(x + 0.5) if x >= 0 else math.ceil(x - 0.5)


def set_formulas(river_fixed: bool, exchange_fixed: bool):
    import importlib
    import rri_torch.hydraulics as hyd
    import rri_torch.exchange as exch
    importlib.reload(hyd)
    importlib.reload(exch)
    river_mod.hq_river = hyd.hq_river if river_fixed else legacy_hq_river
    model_mod.river_slope_exchange = exch.river_slope_exchange if exchange_fixed else legacy_river_slope_exchange


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=2160)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    os.makedirs(RESULTS_DIR, exist_ok=True)

    device = torch.device(args.device)
    dtype = torch.float64

    print("[load] building Solo scenario ...")
    grid, params, outlets, rain = build_solo_scenario(PROJECT_DIR, dtype=dtype)
    dem = read_esri_ascii(os.path.join(PROJECT_DIR, "topo", "adem.txt"))
    rain_seq = build_rain_sequence(rain, grid.nx, grid.ny, dem.xllcorner, dem.yllcorner, dem.cellsize, DT, args.steps, dtype=dtype)

    cepu = outlets["Cepu"]
    idx = cepu
    chain_cpu = [idx]
    while grid.down_riv[idx].item() >= 0:
        idx = grid.down_riv[idx].item()
        chain_cpu.append(idx)
    n_chain = len(chain_cpu)

    grid = grid.to(device)
    for obj, fields in [
        (params.slope, ["ns", "ka", "beta", "da", "dm", "soildepth", "gammaa", "min_wc4latflow"]),
        (params.river, ["width", "depth", "ns_river", "height"]),
    ]:
        for f in fields:
            setattr(obj, f, getattr(obj, f).to(device))
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)
    rain_seq = rain_seq.to(device)

    outlet_names = list(outlets.keys())
    outlet_idx = [outlets[n] for n in outlet_names] + chain_cpu  # gauges first, then the 148-cell chain
    cepu_col = outlet_names.index("Cepu")

    t_ref, q_ref = read_point_hydrograph(os.path.join(PROJECT_DIR, "hydro.txt"))
    stride = int(round(3600.0 / DT))
    n_hours = min(len(t_ref), args.steps // stride)
    sim_idx = np.arange(1, n_hours + 1) * stride - 1
    ref_hourly = q_ref[:n_hours]
    t_hourly = t_ref[:n_hours] / 3600.0

    results = {}
    timings = {}
    for name, cfg in VARIANTS.items():
        set_formulas(cfg["river_fixed"], cfg["exchange_fixed"])
        model = RRIModel(grid, adaptive=True, eps=0.01, ddt_min_slope=1.0, ddt_min_river=0.1, track_qr_avg=True)
        init_state = tuple(x.to(device) for x in model.initial_state(dtype=dtype))
        t0 = time.time()
        out = model.simulate(rain_seq, params, DT, outlet_riv_index=outlet_idx, init_state=init_state)
        if device.type == "cuda":
            torch.cuda.synchronize()
        elapsed = time.time() - t0
        qr_all = out["qr_avg_outlet"].cpu().numpy()
        print(f"[{name}] done in {elapsed:.1f}s, finite={np.isfinite(qr_all).all()}")
        results[name] = qr_all
        timings[name] = elapsed
        np.save(os.path.join(RESULTS_DIR, f"qr_avg_{name}.npy"), qr_all)

    n_gauges = len(outlet_names)

    # --- summary table (Cepu vs hydro.txt) ---
    rows = []
    for name, cfg in VARIANTS.items():
        cepu_series = results[name][:, cepu_col]
        sim_hourly = cepu_series[sim_idx]
        finite = np.isfinite(sim_hourly)
        n = int(np.argmax(~finite)) if not finite.all() else len(finite)
        sim_h, ref_h, t_h = sim_hourly[:n], ref_hourly[:n], t_hourly[:n]
        nse = nash_sutcliffe(ref_h, sim_h)
        rm = rmse(ref_h, sim_h)
        bias = float(np.mean(sim_h - ref_h))
        peak_ref, peak_sim = ref_h.max(), sim_h.max()
        peak_t_ref = t_h[ref_h.argmax()]
        peak_t_sim = t_h[sim_h.argmax()]
        rows.append(dict(
            variant=name, label=cfg["label"], finite_hours=n, nse=nse, rmse=rm, bias=bias,
            peak_ref=peak_ref, peak_sim=peak_sim, peak_t_ref=peak_t_ref, peak_t_sim=peak_t_sim,
            runtime_s=timings[name],
        ))
        print(f"{cfg['label']:32s} NSE={nse:7.4f} RMSE={rm:9.2f} bias={bias:+9.2f} "
              f"peak={peak_sim:9.2f}@{peak_t_sim:.0f}h (ref {peak_ref:.2f}@{peak_t_ref:.0f}h) "
              f"finite_h={n}/{n_hours} runtime={timings[name]:.0f}s")

    csv_path = os.path.join(RESULTS_DIR, "ablation_table.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"[save] {csv_path}")

    md_path = os.path.join(RESULTS_DIR, "ablation_table.md")
    with open(md_path, "w") as f:
        f.write("| Variant | NSE | RMSE (m3/s) | Bias (m3/s) | Peak sim (m3/s) | Peak ref (m3/s) | Runtime (s) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['label']} | {r['nse']:.4f} | {r['rmse']:.1f} | {r['bias']:+.1f} | "
                     f"{r['peak_sim']:.1f} @ {r['peak_t_sim']:.0f}h | {r['peak_ref']:.1f} @ {r['peak_t_ref']:.0f}h | "
                     f"{r['runtime_s']:.0f} |\n")
    print(f"[save] {md_path}")

    # --- figure 1: hydrograph ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t_hourly, ref_hourly, "k-", lw=2.5, label="Official RRI (hydro.txt)", zorder=5)
    for name, cfg in VARIANTS.items():
        cepu_series = results[name][:, cepu_col]
        sim_hourly = cepu_series[sim_idx]
        finite = np.isfinite(sim_hourly)
        n = int(np.argmax(~finite)) if not finite.all() else len(finite)
        ax.plot(t_hourly[:n], sim_hourly[:n], color=COLORS[name], lw=1.6, label=cfg["label"],
                linestyle="-" if name == "both_fixed" else "--")
    ax.set_xlabel("Time [h]")
    ax.set_ylabel(r"Discharge at Cepu [m$^3$/s]")
    ax.set_title("Solo River basin, Cepu gauge: ablation of the two formula fixes")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig_path = os.path.join(RESULTS_DIR, "fig_hydrograph.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    print(f"[save] {fig_path}")

    # --- figure 2: spatial %diff along the chain at t=180h, all 4 variants ---
    step_idx = int(round(180 * 3600 / DT)) - 1  # 0-indexed step for t=180h
    out_dt = args.steps / 96.0
    out_k = None
    for k in range(1, 97):
        if fortran_nint(k * out_dt) - 1 == step_idx:
            out_k = k
            break
    fig2, ax2 = plt.subplots(figsize=(9, 5))
    if out_k is not None and os.path.isdir(FORTRAN_OUT_DIR):
        real_grid = np.loadtxt(os.path.join(FORTRAN_OUT_DIR, f"qr_{out_k:06d}.out"))
        chain_i = [grid.riv_i[k].item() for k in chain_cpu]
        chain_j = [grid.riv_j[k].item() for k in chain_cpu]
        real_vals = np.array([real_grid[i, j] for i, j in zip(chain_i, chain_j)])
        positions = np.arange(n_chain)
        for name, cfg in VARIANTS.items():
            sim_vals = results[name][step_idx, n_gauges:n_gauges + n_chain]
            pct = np.where(real_vals > 1, (sim_vals - real_vals) / np.maximum(np.abs(real_vals), 1.0) * 100, np.nan)
            ax2.plot(positions[:-1], pct[:-1], color=COLORS[name], lw=1.3, label=cfg["label"],
                     linestyle="-" if name == "both_fixed" else "--", alpha=0.85)
        ax2.axhline(0, color="k", lw=0.8)
        ax2.set_xlabel("Position along Cepu -> outlet chain (0 = Cepu)")
        ax2.set_ylabel("% difference vs real Fortran qr_ave")
        ax2.set_title(f"Spatial bias along the river chain at t=180h (outlet cell excluded)")
        ax2.legend(loc="upper right", fontsize=9)
        ax2.grid(alpha=0.3)
        fig2.tight_layout()
        fig2_path = os.path.join(RESULTS_DIR, "fig_spatial_bias.png")
        fig2.savefig(fig2_path, dpi=150)
        print(f"[save] {fig2_path}")
    else:
        print("[skip] fig_spatial_bias.png: real Fortran out/ directory not found or step misaligned")
    plt.close(fig2)

    print("\n[done] ablation complete.")


if __name__ == "__main__":
    main()
