"""The definitive 3-way validation figure for Cepu: official compiled
Fortran simulation, our differentiable rri_torch simulation (both
formula fixes applied -- HANDOFF.md section 7.2), and REAL observed
discharge (`obs/disc_cepu.data`, correctly time-aligned via the
calcHydro output-step-index conversion found in HANDOFF.md section
8.3: `real_hour = column_1 * 3.75`).

This is deliberately just 3 curves, not a multi-model benchmark: the
first two answer "did we port the physics correctly" (they should
nearly overlap), the third answers "does the model match reality"
(expected, honest gap -- this is a validation figure, not a claim that
the port is perfect against the real world). No new simulation is run
here -- reuses `results/qr_avg_both_fixed.npy`, already produced by
`ablation_solo.py`'s full 2160-step run.

Run with:  python examples/final_validation_figure.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from real_data.solo_river import build_solo_scenario, read_point_hydrograph, read_obs_series

PROJECT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "..", "RRI_1_4_2_7_GUI_Beta", "RRI-CUI", "Project", "solo30s")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
DT = 600.0
LASTH, OUTNUM = 360.0, 96.0
STEP_HOURS = LASTH / OUTNUM  # 3.75 -- see HANDOFF.md 8.3


def main():
    # cepu_col within qr_avg_both_fixed.npy's outlet columns -- cheap,
    # no simulation, just re-deriving the same station ordering
    # ablation_solo.py used (outlets dict insertion order).
    _grid, _params, outlets, _rain = build_solo_scenario(PROJECT_DIR)
    outlet_names = list(outlets.keys())
    cepu_col = outlet_names.index("Cepu")

    qr_both_fixed = np.load(os.path.join(RESULTS_DIR, "qr_avg_both_fixed.npy"))
    n_steps = qr_both_fixed.shape[0]
    t_torch = (np.arange(1, n_steps + 1) * DT) / 3600.0  # hours
    q_torch = qr_both_fixed[:, cepu_col]

    t_ref, q_ref = read_point_hydrograph(os.path.join(PROJECT_DIR, "hydro.txt"))
    t_ref = t_ref / 3600.0

    t_obs, q_obs = read_obs_series(os.path.join(PROJECT_DIR, "obs", "disc_cepu.data"), step_hours=STEP_HOURS)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(t_ref, q_ref, "k-", lw=2, label="Official RRI (compiled Fortran)", zorder=3)
    ax.plot(t_torch, q_torch, color="#1f77b4", lw=1.8, ls="--", label="rri_torch (this work, both fixes)", zorder=4)
    ax.plot(t_obs, q_obs, "o", color="#d62728", markersize=7, markeredgecolor="k",
            label="Real observed discharge (obs/disc_cepu.data)", zorder=5)

    ax.set_xlabel("Time [h]")
    ax.set_ylabel(r"Discharge at Cepu [m$^3$/s]")
    ax.set_title("Solo River basin, Cepu gauge: simulation vs. reality")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()

    out_path = os.path.join(RESULTS_DIR, "fig_final_validation.png")
    fig.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    # Quick NSE/RMSE numbers for the caption/text.
    def nse(o, s):
        return 1 - np.sum((o - s) ** 2) / np.sum((o - o.mean()) ** 2)

    def rmse(o, s):
        return float(np.sqrt(np.mean((o - s) ** 2)))

    torch_at_obs = np.interp(t_obs, t_torch, q_torch)
    ref_at_obs = np.interp(t_obs, t_ref, q_ref)
    print(f"rri_torch  vs obs: NSE={nse(q_obs, torch_at_obs):.3f}  RMSE={rmse(q_obs, torch_at_obs):.1f} m3/s")
    print(f"official   vs obs: NSE={nse(q_obs, ref_at_obs):.3f}  RMSE={rmse(q_obs, ref_at_obs):.1f} m3/s")

    torch_at_ref_hours = np.interp(t_ref, t_torch, q_torch)
    print(f"rri_torch  vs official Fortran: NSE={nse(q_ref, torch_at_ref_hours):.4f}  "
          f"RMSE={rmse(q_ref, torch_at_ref_hours):.1f} m3/s")


if __name__ == "__main__":
    main()
