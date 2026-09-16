"""Concrete side-by-side comparison: calibrate the same 4 parameters
against the same synthetic "observed" hydrograph two ways --

  (A) gradient-based: backprop through the differentiable model, Adam.
  (B) gradient-free: treat the model as a black box, Nelder-Mead simplex
      (scipy.optimize), the standard baseline when you can't get gradients.

Both start from the identical wrong initial guess and see the identical
noisy observation. What's compared is how many times the (expensive)
simulator actually has to be run, and how close each gets for that
budget -- this is the concrete version of "gradients let you use way
fewer simulations to calibrate."

Run with:  .venv/Scripts/python.exe examples/compare_gradient_vs_gradientfree.py
"""

from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch
from scipy.optimize import minimize

from rri_torch.model import RRIModel
from examples.synthetic_basin import DTYPE
from examples.calibrate_example import make_scenario

PARAM_NAMES = ["ns_slope", "ka", "ns_river", "ksv"]


def main():
    torch.manual_seed(0)
    grid, true_params, outlet, rain, dt, n_substeps = make_scenario()
    model = RRIModel(grid, n_substeps_slope=n_substeps, n_substeps_river=n_substeps)

    true_values = {
        "ns_slope": true_params.slope.ns[0, 0].item(),
        "ka": true_params.slope.ka[0, 0].item(),
        "ns_river": true_params.river.ns_river.item(),
        "ksv": true_params.ksv[0, 0].item(),
    }
    with torch.no_grad():
        observed = model.simulate(rain, true_params, dt, outlet_riv_index=outlet)["qr_outlet"].clone()
        noise = 0.02 * observed.abs() * torch.randn_like(observed)
        observed = (observed + noise).clamp(min=0.0)

    init_guess = {
        "ns_slope": true_values["ns_slope"] * 0.5,
        "ka": true_values["ka"] * 2.5,
        "ns_river": true_values["ns_river"] * 1.8,
        "ksv": true_values["ksv"] * 0.4,
    }
    x0_log = np.array([np.log(init_guess[k]) for k in PARAM_NAMES])

    def build_params_from_log(log_vec):
        p = true_params
        est = {k: float(np.exp(v)) for k, v in zip(PARAM_NAMES, log_vec)}
        import copy
        p2 = copy.copy(p)
        p2.slope = copy.copy(p.slope)
        p2.river = copy.copy(p.river)
        p2.slope.ns = torch.full_like(p.slope.ns, est["ns_slope"])
        p2.slope.ka = torch.full_like(p.slope.ka, est["ka"])
        p2.river.ns_river = torch.tensor(est["ns_river"], dtype=DTYPE)
        p2.ksv = torch.full_like(p.ksv, est["ksv"])
        return p2, est

    # ---------------------------------------------------------------
    # (B) Gradient-free: Nelder-Mead treats the simulator as a black box.
    # ---------------------------------------------------------------
    nfev_count = [0]

    def loss_np(log_vec):
        nfev_count[0] += 1
        params, _ = build_params_from_log(log_vec)
        with torch.no_grad():
            sim = model.simulate(rain, params, dt, outlet_riv_index=outlet)["qr_outlet"]
        return torch.mean((sim - observed) ** 2).item()

    t0 = time.time()
    result = minimize(
        loss_np, x0_log, method="Nelder-Mead",
        options={"maxfev": 300, "xatol": 1e-6, "fatol": 1e-12, "adaptive": True},
    )
    t_gradfree = time.time() - t0
    _, est_gradfree = build_params_from_log(result.x)

    print("=" * 70)
    print(f"(B) Gradient-free (Nelder-Mead): {nfev_count[0]} simulator calls, "
          f"{t_gradfree:.1f}s, final loss={result.fun:.3e}")
    for k in PARAM_NAMES:
        rel_err = abs(est_gradfree[k] - true_values[k]) / true_values[k]
        print(f"    {k:10s} estimate={est_gradfree[k]:.6g}  true={true_values[k]:.6g}  rel_err={rel_err:.2%}")

    # ---------------------------------------------------------------
    # (A) Gradient-based: Adam using backprop through the simulator.
    # ---------------------------------------------------------------
    log_params = {k: torch.tensor(v, dtype=DTYPE).log().requires_grad_(True) for k, v in init_guess.items()}
    n_iters = 30
    optimizer = torch.optim.Adam(log_params.values(), lr=0.08)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_iters, eta_min=0.005)

    def build_params():
        import copy
        p = copy.copy(true_params)
        p.slope = copy.copy(true_params.slope)
        p.river = copy.copy(true_params.river)
        p.slope.ns = torch.exp(log_params["ns_slope"]).expand_as(true_params.slope.ns)
        p.slope.ka = torch.exp(log_params["ka"]).expand_as(true_params.slope.ka)
        p.river.ns_river = torch.exp(log_params["ns_river"])
        p.ksv = torch.exp(log_params["ksv"]).expand_as(true_params.ksv)
        return p

    best_loss, best_est = float("inf"), None
    t0 = time.time()
    for it in range(n_iters):
        optimizer.zero_grad()
        params = build_params()
        sim = model.simulate(rain, params, dt, outlet_riv_index=outlet)["qr_outlet"]
        loss = torch.mean((sim - observed) ** 2)
        loss.backward()
        optimizer.step()
        scheduler.step()
        if loss.item() < best_loss:
            best_loss = loss.item()
            best_est = {k: torch.exp(v).item() for k, v in log_params.items()}
    t_grad = time.time() - t0

    print("=" * 70)
    print(f"(A) Gradient-based (Adam): {n_iters} simulator calls (each with a free backward pass), "
          f"{t_grad:.1f}s, best loss={best_loss:.3e}")
    for k in PARAM_NAMES:
        rel_err = abs(best_est[k] - true_values[k]) / true_values[k]
        print(f"    {k:10s} estimate={best_est[k]:.6g}  true={true_values[k]:.6g}  rel_err={rel_err:.2%}")

    print("=" * 70)
    print(f"Summary: gradient-free used {nfev_count[0]} simulator calls to reach loss {result.fun:.3e}; "
          f"gradient-based used {n_iters} calls to reach loss {best_loss:.3e}.")


if __name__ == "__main__":
    main()
