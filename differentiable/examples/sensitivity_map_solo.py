"""Item F (plan.md): an adjoint/gradient sensitivity map -- "which grid
cells' Manning's n most affects a downstream discharge" -- computed as
a near-zero-cost byproduct of one backward pass through the real Solo
scenario. This is the concrete illustration of a point unique to
differentiable simulators: a non-differentiable RRI can only tell you
*that* an outlet's flow depends on distant parameters by re-running the
whole simulation once per cell you want to test (18582 extra forward
passes for this domain); one backward pass gives every cell's
sensitivity simultaneously.

Setup mirrors `calibrate_solo.py`'s scoping choices (same file's
docstring has the full reasoning): real Solo geometry, a short strong
synthetic rain pulse (not the real 360h event), and a headwater river
cell as the discharge target (so the "contributing area" that lights up
in the map is meaningful and reachable within a tractable window,
rather than needing 30-40h of river travel time to a far gauge like
Cepu). Uses a longer window than `calibrate_solo.py` (T=80 instead of
T=12) since this is a *single* forward+backward pass, not something
repeated over many optimizer iterations -- affordable at ~10 minutes
(see HANDOFF.md section 8.2's T=80 timing: 614s, 0.60GB peak).

Run with:  python examples/sensitivity_map_solo.py --steps 80 --device cuda
"""
from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch

from rri_torch.model import RRIModel
from real_data.solo_river import build_solo_scenario
from examples.calibrate_solo import find_headwater_cell, to_device

PROJECT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "RRI_1_4_2_7_GUI_Beta", "RRI-CUI", "Project", "solo30s",
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=80)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--rain-mm-hr", type=float, default=60.0)
    args = ap.parse_args()

    device = torch.device(args.device)
    dtype = torch.float64

    grid, params, outlets, _rain_data = build_solo_scenario(PROJECT_DIR)
    grid = grid.to(device)
    params.slope = to_device(params.slope, device)
    params.river = to_device(params.river, device)
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)

    outlet_idx = find_headwater_cell(grid)
    print(f"grid: ny={grid.ny} nx={grid.nx} active={int(grid.domain.sum())}")
    print(f"outlet river-cell index: {outlet_idx} (i={int(grid.riv_i[outlet_idx])}, j={int(grid.riv_j[outlet_idx])})")

    dt = 600.0
    T = args.steps
    rain_rate = args.rain_mm_hr * 1e-3 / 3600.0
    rain_seq = torch.zeros((T, grid.ny, grid.nx), dtype=dtype, device=device)
    pulse_len = max(1, T // 3)
    rain_seq[:pulse_len][:, grid.domain] = rain_rate

    # ns_slope as a full per-cell leaf tensor -- every active cell is its
    # own independent parameter, unlike calibrate_solo.py's single scalar
    # broadcast to the whole domain.
    ns_slope_field = params.slope.ns.clone().requires_grad_(True)
    params.slope.ns = ns_slope_field

    model = RRIModel(grid, n_substeps_slope=40, n_substeps_river=200, checkpoint_substeps=True)
    init_state0 = tuple(x.to(device) for x in model.initial_state(dtype=dtype))

    t0 = time.time()
    out = model.simulate(
        rain_seq, params, dt, outlet_riv_index=outlet_idx, init_state=init_state0, use_checkpointing=True,
    )
    qr = out["qr_outlet"]
    # Total discharged volume through the outlet over the window -- a
    # single scalar target whose sensitivity to every cell's ns_slope we
    # want. (Proportional to sum(qr)*dt; the constant doesn't matter for
    # a *relative* sensitivity map.)
    total_outflow = qr.sum()
    total_outflow.backward()
    print(f"forward+backward: {time.time()-t0:.1f}s, total_outflow={total_outflow.item():.4f} (qr units x steps)")

    grad = ns_slope_field.grad.detach().cpu().numpy()
    domain = grid.domain.cpu().numpy()
    zb = grid.zb.cpu().numpy()
    grad_masked = np.where(domain, grad, np.nan)

    n_nonzero = int((np.abs(grad) > 1e-12) & domain).sum() if False else int(((np.abs(grad) > 1e-12) & domain).sum())
    print(f"cells with nonzero sensitivity: {n_nonzero} / {int(domain.sum())} active cells")
    print(f"|grad| range over nonzero cells: "
          f"[{np.abs(grad[(np.abs(grad)>1e-12)]).min():.3e}, {np.abs(grad).max():.3e}]")

    # ------------------------------------------------------------
    # Figure: sensitivity map (signed, since d(outflow)/d(ns) < 0
    # everywhere it's nonzero -- rougher slope means slower, not more,
    # discharge within a fixed window) with the outlet marked.
    # ------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    ax = axes[0]
    im0 = ax.imshow(np.where(domain, zb, np.nan), cmap="terrain")
    oi, oj = int(grid.riv_i[outlet_idx]), int(grid.riv_j[outlet_idx])
    ax.plot(oj, oi, "r*", markersize=16, markeredgecolor="k", label="outlet (headwater cell)")
    ax.set_title("Elevation (context)")
    ax.legend(loc="upper right", fontsize=8)
    plt.colorbar(im0, ax=ax, label="m", fraction=0.046)

    ax = axes[1]
    vmax = np.nanpercentile(np.abs(grad_masked), 99)
    im1 = ax.imshow(grad_masked, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.plot(oj, oi, "k*", markersize=16, markeredgecolor="w")
    ax.set_title(f"d(total outflow)/d(ns_slope) per cell\n(T={T} steps, one backward pass)")
    plt.colorbar(im1, ax=ax, label="sensitivity [m^3/s per unit ns]", fraction=0.046)

    for ax in axes:
        ax.set_xlabel("column")
        ax.set_ylabel("row")

    fig.suptitle("Solo River basin: adjoint sensitivity of headwater discharge to spatial Manning's n")
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "fig_sensitivity_map.png"), dpi=150)
    print(f"saved {out_dir}/fig_sensitivity_map.png")

    np.savez(os.path.join(out_dir, "sensitivity_map_solo.npz"), grad=grad, domain=domain, zb=zb,
             outlet_i=oi, outlet_j=oj, steps=T)
    print(f"saved {out_dir}/sensitivity_map_solo.npz")


if __name__ == "__main__":
    main()
