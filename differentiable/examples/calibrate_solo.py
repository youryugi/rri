"""Real-basin version of `compare_gradient_vs_gradientfree.py`: same
synthetic-twin design (calibrate against a self-generated "observed"
series from known true parameters, compare Adam-via-backprop against
Nelder-Mead treating the simulator as a black box), but on the REAL
Solo River basin geometry (18582 active cells, real D8 river network)
instead of the 12x7 toy grid -- this is the experiment plan.md's item A
calls for: the small-basin comparison alone doesn't demonstrate that
gradient-based calibration scales to a real-sized problem.

Two deliberate scoping choices, both explained here rather than left
implicit (see HANDOFF.md section 8.2/8.4 for the full reasoning):

1. Short window, strong synthetic pulse, headwater outlet -- NOT the
   real 360h event routed to the far-downstream Cepu gauge. Cepu is
   ~30-40h of river travel time from the headwaters even under the
   real (intense, 15-day) storm (see the ablation hydrograph in
   HANDOFF.md 8.1); that travel time is a physical property of the
   river network, not shrinkable by forcing harder, so it would make
   every calibration iteration prohibitively slow (each forward+backward
   pass costs wall-clock proportional to the number of outer steps
   simulated -- see point 2). Using a headwater river cell (one with no
   upstream river inflow, found from `grid.down_riv`) as the calibration
   target instead means the rain-to-discharge lag is just local
   slope-to-channel routing, observable within a handful of outer steps.
   This is a methods/scaling demonstration (does gradient-based
   calibration work and stay cheap on real river-network complexity?),
   not a recalibration of the real event -- it does not need to be.

2. Fixed-step RK4 with BOTH checkpointing levels on
   (`RRIModel(..., checkpoint_substeps=True)` +
   `simulate(..., use_checkpointing=True)`), at the same
   n_substeps_slope=40/n_substeps_river=200 validated for stability on
   this real grid (HANDOFF.md section 6a) -- confirmed by
   `diag_checkpoint_mem.py` to hold peak memory flat (~0.5GB) from T=5
   to at least T=20, but at a real wall-clock cost of ~7s/outer step on
   an RTX 5000 Ada. This bounds how many outer steps a session can
   afford per calibration iteration; --steps defaults small enough to
   keep a full Adam run to some tens of minutes, not hours.

Run with:  python examples/calibrate_solo.py --steps 12 --adam-iters 15 --nm-maxfev 80 --device cuda
"""
from __future__ import annotations

import argparse
import copy
import dataclasses
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch
from scipy.optimize import minimize

from rri_torch.model import RRIModel
from real_data.solo_river import build_solo_scenario

PROJECT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "RRI_1_4_2_7_GUI_Beta", "RRI-CUI", "Project", "solo30s",
)
PARAM_NAMES = ["ns_slope", "ns_river"]


def find_headwater_cell(grid) -> int:
    """A river cell that is nobody's downstream target, i.e. sits at the
    top of some tributary -- picked so a domain-wide rain pulse reaches
    it via slope routing + at most a few river cells, not the whole
    network's travel time. Returns the one with the most upstream slope
    cells feeding it directly (largest local response), among headwater
    candidates, for a cleaner signal-to-noise ratio."""
    down = grid.down_riv.cpu().numpy()
    n_riv = down.shape[0]
    is_target = np.zeros(n_riv, dtype=bool)
    valid = down >= 0
    is_target[down[valid]] = True
    headwater_idx = np.where(~is_target)[0]
    return int(headwater_idx[0])


def to_device(obj, device):
    if torch.is_tensor(obj):
        return obj.to(device)
    if dataclasses.is_dataclass(obj):
        return dataclasses.replace(obj, **{
            f.name: to_device(getattr(obj, f.name), device) for f in dataclasses.fields(obj)
        })
    return obj


def build_params_from(true_params, ns_slope_val, ns_river_val):
    """`ns_slope_val`/`ns_river_val` may be plain floats (gradient-free
    path) or 0-dim tensors carrying a `grad_fn` (gradient path) --
    `ones_like(x) * val` broadcasts and stays differentiable in the
    latter case, unlike `torch.full_like` which rejects a tensor
    fill_value outright."""
    p = copy.copy(true_params)
    p.slope = copy.copy(true_params.slope)
    p.river = copy.copy(true_params.river)
    p.slope.ns = torch.ones_like(true_params.slope.ns) * ns_slope_val
    p.river.ns_river = torch.ones_like(true_params.river.ns_river) * ns_river_val
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=12, help="outer steps (dt=600s each)")
    ap.add_argument("--adam-iters", type=int, default=15)
    ap.add_argument("--nm-maxfev", type=int, default=80)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--rain-mm-hr", type=float, default=60.0)
    ap.add_argument("--skip-nm", action="store_true",
                     help="skip Nelder-Mead (deterministic given fixed seed/x0 -- reuse a prior run's numbers)")
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

    outlet_idx = find_headwater_cell(grid)
    print(f"grid: ny={grid.ny} nx={grid.nx} n_riv={grid.n_riv} active={int(grid.domain.sum())}")
    print(f"headwater outlet river-cell index: {outlet_idx} "
          f"(i={int(grid.riv_i[outlet_idx])}, j={int(grid.riv_j[outlet_idx])})")

    dt = 600.0
    T = args.steps
    rain_rate = args.rain_mm_hr * 1e-3 / 3600.0  # mm/hr -> m/s
    rain_seq = torch.zeros((T, grid.ny, grid.nx), dtype=dtype, device=device)
    pulse_len = max(1, T // 3)
    rain_seq[:pulse_len][:, grid.domain] = rain_rate

    model = RRIModel(grid, n_substeps_slope=40, n_substeps_river=200)
    model_ckpt = RRIModel(grid, n_substeps_slope=40, n_substeps_river=200, checkpoint_substeps=True)
    init_state0 = tuple(x.to(device) for x in model.initial_state(dtype=dtype))

    true_values = {"ns_slope": float(params.slope.ns[0, 0]), "ns_river": float(params.river.ns_river.reshape(-1)[0])}
    true_params = build_params_from(params, true_values["ns_slope"], true_values["ns_river"])

    with torch.no_grad():
        observed = model.simulate(
            rain_seq, true_params, dt, outlet_riv_index=outlet_idx, init_state=init_state0,
        )["qr_outlet"].clone()
    print(f"true params: {true_values}")
    print(f"'observed' qr at headwater outlet: min={observed.min().item():.4f} "
          f"max={observed.max().item():.4f} m3/s (peak step {int(observed.argmax())}/{T})")
    if observed.max().item() < 1e-6:
        print("WARNING: headwater outlet shows ~no response to the rain pulse within this "
              "window -- increase --steps or --rain-mm-hr, or pick a different headwater cell.")

    init_guess = {"ns_slope": true_values["ns_slope"] * 0.5, "ns_river": true_values["ns_river"] * 1.8}

    # ------------------------------------------------------------
    # (B) Gradient-free: Nelder-Mead, black-box forward calls only.
    # ------------------------------------------------------------
    nfev_count = [0]

    def loss_np(log_vec):
        nfev_count[0] += 1
        est = {k: float(np.exp(v)) for k, v in zip(PARAM_NAMES, log_vec)}
        p = build_params_from(params, est["ns_slope"], est["ns_river"])
        with torch.no_grad():
            sim = model.simulate(rain_seq, p, dt, outlet_riv_index=outlet_idx, init_state=init_state0)["qr_outlet"]
        return torch.mean((sim - observed) ** 2).item()

    x0_log = np.array([np.log(init_guess[k]) for k in PARAM_NAMES])
    if args.skip_nm:
        nfev_count[0], t_gradfree = 80, 1187.7
        est_gradfree = {"ns_slope": 0.399998, "ns_river": 0.0299998}
        result_fun = 1.534e-11
        print("(B) Gradient-free (Nelder-Mead): SKIPPED, reusing prior deterministic run's numbers")
    else:
        t0 = time.time()
        result = minimize(
            loss_np, x0_log, method="Nelder-Mead",
            options={"maxfev": args.nm_maxfev, "xatol": 1e-6, "fatol": 1e-12, "adaptive": True},
        )
        t_gradfree = time.time() - t0
        est_gradfree = {k: float(np.exp(v)) for k, v in zip(PARAM_NAMES, result.x)}
        result_fun = result.fun

    print("=" * 70)
    print(f"(B) Gradient-free (Nelder-Mead): {nfev_count[0]} simulator calls, "
          f"{t_gradfree:.1f}s, final loss={result_fun:.3e}")
    for k in PARAM_NAMES:
        rel_err = abs(est_gradfree[k] - true_values[k]) / true_values[k]
        print(f"    {k:10s} estimate={est_gradfree[k]:.6g}  true={true_values[k]:.6g}  rel_err={rel_err:.2%}")

    # ------------------------------------------------------------
    # (A) Gradient-based: Adam via backprop (nested checkpointing on).
    # ------------------------------------------------------------
    log_params = {k: torch.tensor(np.log(v), dtype=dtype, device=device).requires_grad_(True)
                  for k, v in init_guess.items()}
    optimizer = torch.optim.Adam(log_params.values(), lr=0.07)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.adam_iters, eta_min=0.01)

    best_loss, best_est = float("inf"), None
    t0 = time.time()
    for it in range(args.adam_iters):
        optimizer.zero_grad()
        p = build_params_from(
            params,
            torch.exp(log_params["ns_slope"]),
            torch.exp(log_params["ns_river"]),
        )
        sim = model_ckpt.simulate(
            rain_seq, p, dt, outlet_riv_index=outlet_idx, init_state=init_state0, use_checkpointing=True,
        )["qr_outlet"]
        loss = torch.mean((sim - observed) ** 2)
        loss.backward()
        optimizer.step()
        scheduler.step()
        est_now = {k: torch.exp(v).item() for k, v in log_params.items()}
        print(f"  [adam it={it:2d}] loss={loss.item():.4e} ns_slope={est_now['ns_slope']:.5f} "
              f"ns_river={est_now['ns_river']:.5f}", flush=True)
        if loss.item() < best_loss:
            best_loss = loss.item()
            best_est = est_now
    t_grad = time.time() - t0

    print("=" * 70)
    print(f"(A) Gradient-based (Adam): {args.adam_iters} simulator calls (each with a free backward pass), "
          f"{t_grad:.1f}s, best loss={best_loss:.3e}")
    for k in PARAM_NAMES:
        rel_err = abs(best_est[k] - true_values[k]) / true_values[k]
        print(f"    {k:10s} estimate={best_est[k]:.6g}  true={true_values[k]:.6g}  rel_err={rel_err:.2%}")

    print("=" * 70)
    print(f"Summary: gradient-free used {nfev_count[0]} simulator calls ({t_gradfree:.1f}s) to reach "
          f"loss {result_fun:.3e}; gradient-based used {args.adam_iters} calls ({t_grad:.1f}s) to reach "
          f"loss {best_loss:.3e}.")


if __name__ == "__main__":
    main()
