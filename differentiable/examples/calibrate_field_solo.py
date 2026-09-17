"""Item E (plan.md): calibrate a full spatially-distributed Manning's n
FIELD (one free parameter per active grid cell, ~18582 of them here)
instead of `calibrate_solo.py`'s two scalars -- this is the experiment
that makes the "why differentiable" argument concrete: Nelder-Mead-style
derivative-free search is not just slower but outright infeasible at
this dimensionality (a simplex in 18582 dimensions needs >18582 points
just to be non-degenerate), while one backward pass gives every cell's
gradient at once, at the same cost as `calibrate_solo.py`'s 2-parameter
version.

Design notes (see HANDOFF.md section 8.5 for full writeup):

- **True field**: `ns_true = 0.2 + 0.5 * normalized_elevation` -- a
  physically-plausible story (higher/steeper ground tends to be
  rougher: more vegetation, less channelization) and, more importantly
  for demonstrating recovery, a *smooth, structured* field rather than
  per-cell noise, since with a handful of point observations this
  inverse problem is severely underdetermined without some structural
  assumption tying nearby cells together (see the smoothness
  regularizer below).

- **Multiple observation points**: a single scalar time series (one
  outlet) cannot possibly constrain 18582 independent parameters --
  this uses `N_OUTLETS` headwater river cells (see
  `calibrate_solo.py.find_headwater_cell`, extended here to return all
  of them) spread across the basin, so different sub-catchments each
  get their own constraint. This is still far short of one observation
  per parameter -- see the honest identifiability discussion in the
  results writeup: recovery quality is reported *separately* for cells
  within reach of some observation window vs. cells that are not (the
  latter are, correctly, unconstrained and stay near the initial guess
  -- exactly the behavior a real distributed calibration would show,
  and worth reporting rather than hiding).

- **Smoothness regularization**: an L2 penalty on `ns` differences
  between orthogonally-adjacent active cells. Without it the inverse
  problem is so underdetermined that Adam can drive isolated,
  disconnected cells to extreme values while barely moving the loss
  (classic ill-posed-inverse-problem overfitting); the penalty encodes
  the (true, here) prior that roughness varies smoothly, not
  cell-to-cell.

Run with:  python examples/calibrate_field_solo.py --steps 30 --iters 15 --device cuda
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
from examples.calibrate_solo import to_device

PROJECT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "RRI_1_4_2_7_GUI_Beta", "RRI-CUI", "Project", "solo30s",
)

# 8 headwater river cells spread across the basin (row/col range
# 33-148 / 42-277 out of the full 69 headwater cells -- see
# HANDOFF.md 8.5 for how these were picked), used as independent
# multi-point observation targets.
OUTLET_INDICES = [0, 190, 902, 930, 567, 298, 630, 813]


def find_all_headwater_cells(grid) -> np.ndarray:
    down = grid.down_riv.cpu().numpy()
    is_target = np.zeros(down.shape[0], dtype=bool)
    valid = down >= 0
    is_target[down[valid]] = True
    return np.where(~is_target)[0]


def smoothness_penalty(field: torch.Tensor, domain: torch.Tensor) -> torch.Tensor:
    """Sum of squared differences between orthogonally-adjacent active
    cells (both endpoints must be in-domain, so this never regularizes
    across the basin boundary into meaningless exterior values)."""
    both_h = domain[:, :-1] & domain[:, 1:]
    both_v = domain[:-1, :] & domain[1:, :]
    dh = (field[:, :-1] - field[:, 1:])[both_h]
    dv = (field[:-1, :] - field[1:, :])[both_v]
    return (dh ** 2).sum() + (dv ** 2).sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=30)
    ap.add_argument("--iters", type=int, default=15)
    ap.add_argument("--lr", type=float, default=0.05)
    ap.add_argument("--reg-weight", type=float, default=2.0)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--rain-mm-hr", type=float, default=60.0)
    args = ap.parse_args()

    device = torch.device(args.device)
    dtype = torch.float64
    torch.manual_seed(0)

    grid, params, outlets, _rain_data = build_solo_scenario(PROJECT_DIR)
    grid = grid.to(device)
    params.slope = to_device(params.slope, device)
    params.river = to_device(params.river, device)
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)

    domain = grid.domain
    zb = grid.zb
    zb_active = zb[domain]
    zb_norm = torch.zeros_like(zb)
    zb_norm[domain] = (zb_active - zb_active.min()) / (zb_active.max() - zb_active.min())
    ns_true_field = torch.where(domain, 0.2 + 0.5 * zb_norm, torch.zeros_like(zb))

    dt = 600.0
    T = args.steps
    rain_rate = args.rain_mm_hr * 1e-3 / 3600.0
    rain_seq = torch.zeros((T, grid.ny, grid.nx), dtype=dtype, device=device)
    pulse_len = max(1, T // 3)
    rain_seq[:pulse_len][:, domain] = rain_rate

    model_plain = RRIModel(grid, n_substeps_slope=40, n_substeps_river=200)
    model_ckpt = RRIModel(grid, n_substeps_slope=40, n_substeps_river=200, checkpoint_substeps=True)
    init_state0 = tuple(x.to(device) for x in model_plain.initial_state(dtype=dtype))

    outlet_idx = torch.tensor(OUTLET_INDICES, device=device)
    print(f"grid: ny={grid.ny} nx={grid.nx} active={int(domain.sum())}")
    print(f"{len(OUTLET_INDICES)} observation outlets at "
          f"{[(int(grid.riv_i[i]), int(grid.riv_j[i])) for i in OUTLET_INDICES]}")

    params_true = params
    params_true.slope.ns = ns_true_field
    with torch.no_grad():
        observed = model_plain.simulate(
            rain_seq, params_true, dt, outlet_riv_index=outlet_idx, init_state=init_state0,
        )["qr_outlet"].clone()
    print(f"'observed' peak per outlet: {observed.max(dim=0).values.tolist()}")

    # Which cells actually influence at least one observation within
    # this window -- the honest identifiability check. Cheap to get:
    # reuse the sensitivity-map trick (one backward pass of the *true*
    # field's outflow sum w.r.t. ns) rather than re-deriving reachability
    # from the D8 tree by hand.
    ns_probe = ns_true_field.clone().requires_grad_(True)
    params_probe = params
    params_probe.slope.ns = ns_probe
    probe_out = model_ckpt.simulate(
        rain_seq, params_probe, dt, outlet_riv_index=outlet_idx, init_state=init_state0, use_checkpointing=True,
    )["qr_outlet"].sum()
    probe_out.backward()
    reachable = (ns_probe.grad.detach().abs() > 1e-12) & domain
    print(f"cells reachable from at least one outlet within T={T}: {int(reachable.sum())} / {int(domain.sum())}")

    # ------------------------------------------------------------
    # Calibrate the full field from a flat, wrong initial guess.
    # ------------------------------------------------------------
    ns_init = torch.full_like(zb, 0.4)
    log_field = torch.log(ns_init.clamp(min=1e-3)).clone().requires_grad_(True)
    optimizer = torch.optim.Adam([log_field], lr=args.lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.iters, eta_min=args.lr * 0.1)

    loss_history = []
    t0 = time.time()
    for it in range(args.iters):
        optimizer.zero_grad()
        ns_field = torch.exp(log_field)
        params.slope.ns = ns_field
        sim = model_ckpt.simulate(
            rain_seq, params, dt, outlet_riv_index=outlet_idx, init_state=init_state0, use_checkpointing=True,
        )["qr_outlet"]
        data_loss = torch.mean((sim - observed) ** 2)
        reg = smoothness_penalty(ns_field, domain) / max(1, int(domain.sum()))
        loss = data_loss + args.reg_weight * reg
        loss.backward()
        optimizer.step()
        scheduler.step()
        loss_history.append((data_loss.item(), reg.item()))
        print(f"  [it={it:2d}] data_loss={data_loss.item():.4e} reg={reg.item():.4e} "
              f"total={loss.item():.4e}", flush=True)
    t_total = time.time() - t0

    with torch.no_grad():
        ns_recovered = torch.exp(log_field).detach()

    err = (ns_recovered - ns_true_field).cpu().numpy()
    reach_np = reachable.cpu().numpy()
    domain_np = domain.cpu().numpy()
    rmse_reachable = float(np.sqrt(np.mean(err[reach_np] ** 2)))
    rmse_unreachable = float(np.sqrt(np.mean(err[domain_np & ~reach_np] ** 2))) if (domain_np & ~reach_np).any() else float("nan")
    print("=" * 70)
    print(f"Calibration done: {args.iters} iterations, {t_total:.1f}s")
    print(f"RMSE(ns) over cells reachable from an outlet: {rmse_reachable:.4f}")
    print(f"RMSE(ns) over cells NOT reachable (should stay near init=0.4, "
          f"true ranges 0.2-0.7): {rmse_unreachable:.4f}")

    # ------------------------------------------------------------
    # Figure: true field, recovered field, error map.
    # ------------------------------------------------------------
    zb_np = zb.cpu().numpy()
    true_np = ns_true_field.cpu().numpy()
    rec_np = ns_recovered.cpu().numpy()

    fig, axes = plt.subplots(1, 4, figsize=(20, 5.5))
    vmin, vmax = 0.2, 0.7

    ax = axes[0]
    im = ax.imshow(np.where(domain_np, zb_np, np.nan), cmap="terrain")
    plt.colorbar(im, ax=ax, label="m", fraction=0.046)
    ax.set_title("Elevation (context)")

    ax = axes[1]
    im = ax.imshow(np.where(domain_np, true_np, np.nan), cmap="viridis", vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=ax, label="ns_slope", fraction=0.046)
    ax.set_title("True ns_slope field")

    ax = axes[2]
    im = ax.imshow(np.where(domain_np, rec_np, np.nan), cmap="viridis", vmin=vmin, vmax=vmax)
    for oi in OUTLET_INDICES:
        ax.plot(int(grid.riv_j[oi]), int(grid.riv_i[oi]), "r*", markersize=10, markeredgecolor="k")
    ax.set_title(f"Recovered field ({args.iters} Adam iters)\nred stars = observation outlets")
    plt.colorbar(im, ax=ax, label="ns_slope", fraction=0.046)

    ax = axes[3]
    err_masked = np.where(domain_np, err, np.nan)
    vabs = np.nanpercentile(np.abs(err_masked), 95)
    im = ax.imshow(err_masked, cmap="RdBu_r", vmin=-vabs, vmax=vabs)
    ax.contourf(reach_np.astype(float), levels=[0.5, 1.5], colors="none", hatches=["..."], alpha=0)
    plt.colorbar(im, ax=ax, label="recovered - true", fraction=0.046)
    ax.set_title("Error (recovered - true)")

    for ax in axes:
        ax.set_xlabel("column")
        ax.set_ylabel("row")

    fig.suptitle("Solo River basin: spatially-distributed Manning's n calibration from 8 headwater outlets")
    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, "fig_field_calibration.png"), dpi=150)
    print(f"saved {out_dir}/fig_field_calibration.png")

    np.savez(os.path.join(out_dir, "field_calibration_solo.npz"), true=true_np, recovered=rec_np,
             domain=domain_np, reachable=reach_np, zb=zb_np, loss_history=np.array(loss_history))
    print(f"saved {out_dir}/field_calibration_solo.npz")


if __name__ == "__main__":
    main()
