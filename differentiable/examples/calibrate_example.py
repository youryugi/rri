"""Gradient-based calibration demo.

Generates a synthetic "observed" outlet hydrograph from a catchment with
known parameters, then recovers those parameters from a different, wrong
starting guess purely by backpropagating an MSE loss through the whole
routing model and taking Adam steps -- no forward-difference finite
differencing, no surrogate, just autograd through the physics.

This is the point of making RRI differentiable at all: this exact loop is
what you'd run against a real observed hydrograph to calibrate ns_slope,
ka, ns_river, ksv, etc. for an actual catchment.

Run with:  .venv/Scripts/python.exe examples/calibrate_example.py
"""

from __future__ import annotations

import copy
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch

from rri_torch.model import RRIModel
from examples.synthetic_basin import build_synthetic_basin, make_rain_pulse, DTYPE


def make_scenario():
    grid, params, outlet = build_synthetic_basin()
    # A shallow soil layer + a heavy storm exercises both the subsurface
    # (ka, beta) and Manning surface-flow (ns_slope) branches, and the
    # river branch (ns_river) always participates -- so every parameter
    # calibrated below actually influences the outlet hydrograph.
    soildepth = torch.full_like(params.slope.soildepth, 0.05)
    da = soildepth * params.slope.gammaa
    dm = 0.3 * da
    params.slope.soildepth = soildepth
    params.slope.da = da
    params.slope.dm = dm

    T = 16
    rain = make_rain_pulse(T, grid.ny, grid.nx, peak_mm_per_hr=80.0, start=2, duration=6)
    dt = 300.0
    n_substeps = 10
    return grid, params, outlet, rain, dt, n_substeps


def main():
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

    # --- calibration setup: log-parametrize so Adam steps stay positive
    # and comparable in scale across parameters spanning several orders
    # of magnitude (ka ~ 1e-5, ns_slope ~ 1e-1, ...).
    init_guess = {
        "ns_slope": true_values["ns_slope"] * 0.5,
        "ka": true_values["ka"] * 2.5,
        "ns_river": true_values["ns_river"] * 1.8,
        "ksv": true_values["ksv"] * 0.4,
    }
    log_params = {k: torch.tensor(v, dtype=DTYPE).log().requires_grad_(True) for k, v in init_guess.items()}

    n_iters = 30
    optimizer = torch.optim.Adam(log_params.values(), lr=0.08)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_iters, eta_min=0.005)

    def build_params():
        p = copy.copy(true_params)
        p.slope = copy.copy(true_params.slope)
        p.river = copy.copy(true_params.river)
        p.slope.ns = torch.exp(log_params["ns_slope"]).expand_as(true_params.slope.ns)
        p.slope.ka = torch.exp(log_params["ka"]).expand_as(true_params.slope.ka)
        p.river.ns_river = torch.exp(log_params["ns_river"])
        p.ksv = torch.exp(log_params["ksv"]).expand_as(true_params.ksv)
        return p

    best_loss = float("inf")
    best_est = None
    t0 = time.time()
    for it in range(n_iters):
        optimizer.zero_grad()
        params = build_params()
        sim = model.simulate(rain, params, dt, outlet_riv_index=outlet)["qr_outlet"]
        loss = torch.mean((sim - observed) ** 2)
        loss.backward()
        optimizer.step()
        scheduler.step()

        loss_v = loss.item()
        if loss_v < best_loss:
            best_loss = loss_v
            best_est = {k: torch.exp(v).item() for k, v in log_params.items()}

        if it % 5 == 0 or it == n_iters - 1:
            est = {k: torch.exp(v).item() for k, v in log_params.items()}
            print(
                f"iter {it:3d}  loss={loss_v:.6e}  "
                + "  ".join(f"{k}={est[k]:.4g}(true {true_values[k]:.4g})" for k in est)
            )

    print(f"\nDone in {time.time() - t0:.1f}s. Best-loss estimates ({best_loss:.3e}) vs. truth:")
    for k, est in best_est.items():
        true_v = true_values[k]
        print(f"  {k:10s} estimate={est:.6g}  true={true_v:.6g}  rel_err={abs(est - true_v) / true_v:.2%}")


if __name__ == "__main__":
    main()
