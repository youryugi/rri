"""Compares rri_torch's Cepu discharge hydrograph (from a saved
solo_validation.py --save run) against RRI's own pre-computed reference,
`solo30s/hydro.txt` -- a point time series written directly by RRI.f90
itself during ICHARM's official run (hourly, NOT derived from the
qr_*.out spatial dumps; see real_data/solo_river.read_point_hydrograph).

This does not require re-running the Fortran binary: `hydro.txt` already
ships in the project as ICHARM's own reference output.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch

from real_data.solo_river import read_point_hydrograph

PROJECT_DIR = "/home/yang/github/rri/RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s"
DT = 600.0


def nash_sutcliffe(obs, sim):
    return 1.0 - np.sum((obs - sim) ** 2) / np.sum((obs - obs.mean()) ** 2)


def rmse(obs, sim):
    return float(np.sqrt(np.mean((obs - sim) ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="path to a solo_validation.py --save .pt file")
    ap.add_argument("--instantaneous", action="store_true",
                     help="compare against the instantaneous end-of-step qr even if qr_avg was saved "
                          "(default: prefer qr_avg, since that's what hydro.txt's qr_ave actually is)")
    args = ap.parse_args()

    d = torch.load(args.run, weights_only=False)
    names = d["names"]
    use_avg = ("qr_avg" in d) and not args.instantaneous
    qr = (d["qr_avg"] if use_avg else d["qr"]).numpy()  # (T, n_gauges)
    print(f"[sim] comparing against {'time-averaged qr_avg (RRI qr_ave-equivalent)' if use_avg else 'instantaneous end-of-step qr'}")
    if "Cepu" not in names:
        print("Cepu not among saved gauges:", names)
        return
    cepu_idx = names.index("Cepu")
    qr_cepu = qr[:, cepu_idx]
    n_steps = qr_cepu.shape[0]

    t_ref, q_ref = read_point_hydrograph(os.path.join(PROJECT_DIR, "hydro.txt"))
    print(f"[ref] hydro.txt: {len(t_ref)} hourly points, t in [{t_ref[0]:.0f}, {t_ref[-1]:.0f}]s")
    print(f"[sim] rri_torch: {n_steps} steps @ dt={DT}s, t in [{DT:.0f}, {n_steps*DT:.0f}]s")

    # rri_torch's step k (0-indexed) covers the window [k*dt, (k+1)*dt] and
    # is indexed by its end time (k+1)*dt -- qr[k] is instantaneous at that
    # instant, qr_avg[k] is the ddt-weighted time-average *over* that same
    # window (matching what RRI.f90 actually resets and re-accumulates
    # `qr_ave` over every outer step, not over the full reporting hour --
    # see HANDOFF.md section 6a). hydro.txt is hourly starting at t=3600s,
    # so hourly point m (1-indexed) <-> sim step index (m*3600/dt - 1).
    stride = int(round(3600.0 / DT))
    n_hours = min(len(t_ref), n_steps // stride)
    sim_idx = np.arange(1, n_hours + 1) * stride - 1
    sim_hourly = qr_cepu[sim_idx]
    ref_hourly = q_ref[:n_hours]
    t_hourly = t_ref[:n_hours]

    finite = np.isfinite(sim_hourly)
    print(f"[compare] {n_hours} hourly points compared, {finite.sum()} finite in sim")
    if finite.sum() < n_hours:
        first_bad = int(np.argmax(~finite))
        print(f"[compare] first non-finite at hour index {first_bad} (t={t_hourly[first_bad]:.0f}s = "
              f"{t_hourly[first_bad]/3600:.1f}h)")
        n_hours = first_bad
        sim_hourly = sim_hourly[:n_hours]
        ref_hourly = ref_hourly[:n_hours]
        t_hourly = t_hourly[:n_hours]

    nse = nash_sutcliffe(ref_hourly, sim_hourly)
    rmse_v = rmse(ref_hourly, sim_hourly)
    bias = float(np.mean(sim_hourly - ref_hourly))
    peak_ref, peak_sim = ref_hourly.max(), sim_hourly.max()
    peak_t_ref = t_hourly[ref_hourly.argmax()] / 3600.0
    peak_t_sim = t_hourly[sim_hourly.argmax()] / 3600.0

    print(f"\n=== Cepu gauge: rri_torch vs RRI Fortran (hydro.txt) over {n_hours}h ===")
    print(f"NSE   = {nse:.4f}")
    print(f"RMSE  = {rmse_v:.3f} m^3/s")
    print(f"bias  = {bias:+.3f} m^3/s (mean sim - ref)")
    print(f"peak: ref={peak_ref:.2f} m^3/s @ {peak_t_ref:.1f}h   sim={peak_sim:.2f} m^3/s @ {peak_t_sim:.1f}h")

    print(f"\n{'t[h]':>8s} {'ref[m3/s]':>12s} {'sim[m3/s]':>12s} {'diff':>10s}")
    step = max(1, n_hours // 24)
    for i in range(0, n_hours, step):
        print(f"{t_hourly[i]/3600:8.1f} {ref_hourly[i]:12.3f} {sim_hourly[i]:12.3f} {sim_hourly[i]-ref_hourly[i]:10.3f}")


if __name__ == "__main__":
    main()
