"""Plain-script sanity checks (no pytest dependency):

1. The model runs forward on the synthetic basin without NaNs.
2. Water mass is (approximately) conserved: rainfall in ~= outlet
   discharge out + storage change + infiltration + drained volume.
3. Gradients of the outlet hydrograph w.r.t. key learnable parameters
   (slope Manning's n, lateral conductivity, river Manning's n, Green-Ampt
   Ks) are finite and non-zero.

Run with:  .venv/Scripts/python.exe tests/test_forward_and_grad.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch

from rri_torch.model import RRIModel
from examples.synthetic_basin import build_synthetic_basin, make_rain_pulse, DTYPE


def test_forward_runs():
    grid, params, outlet = build_synthetic_basin()
    model = RRIModel(grid, n_substeps_slope=4, n_substeps_river=4)
    T = 40
    rain = make_rain_pulse(T, grid.ny, grid.nx)
    dt = 300.0  # 5 minutes

    out = model.simulate(rain, params, dt, outlet_riv_index=outlet)
    qr = out["qr_outlet"]

    assert torch.isfinite(qr).all(), "non-finite discharge in hydrograph"
    assert torch.isfinite(out["hs_final"]).all()
    assert torch.isfinite(out["hr_final"]).all()
    assert qr.max().item() > 0.0, "expected the pulse to produce outlet discharge"
    print(f"[forward] peak outlet discharge = {qr.max().item():.4f} m^3/s at step {int(qr.argmax())}/{T}")
    return grid, params, model, rain, dt, outlet, out


def test_mass_balance(grid, params, model, rain, dt, outlet):
    # Re-run keeping track of everything needed for a mass-balance check.
    hs, gampt_ff, hr = model.initial_state(dtype=DTYPE)
    total_rain_in = 0.0
    total_infilt = 0.0
    total_drained = 0.0
    T = rain.shape[0]
    for t in range(T):
        active_cells = grid.domain.sum().item()
        total_rain_in += rain[t][grid.domain].sum().item() * grid.area * dt
        hs, gampt_ff, hr, diag = model.step(hs, gampt_ff, hr, params, rain[t], None, dt)
        total_infilt += diag.infilt_rate[grid.domain].sum().item() * grid.area * dt
        total_drained += diag.drained_volume.item()

    storage_slope = (hs[grid.domain] * grid.area).sum().item()
    storage_infilt = (gampt_ff[grid.domain] * grid.area).sum().item()
    channel_area = (params.river.width * grid.len_riv)
    storage_river = (hr * channel_area).sum().item()

    balance_in = total_rain_in
    balance_out = storage_slope + storage_infilt + storage_river + total_drained
    # total_infilt is not an independent loss term: it just moves water
    # from `hs` into `gampt_ff`, both already counted in storage above.
    rel_err = abs(balance_in - balance_out) / max(balance_in, 1e-9)
    print(f"[mass balance] rain_in={balance_in:.3e} m^3, "
          f"stored+drained={balance_out:.3e} m^3, rel_err={rel_err:.3%}")
    assert rel_err < 0.05, f"mass balance relative error too large: {rel_err:.3%}"


def test_gradients():
    grid, params, outlet = build_synthetic_basin()
    # This scenario (thin soil layer, heavy storm, see below) is numerically
    # stiffer than the smoke test above -- the reference Fortran handles
    # that by shrinking its adaptive RK step; we use a fixed step, so we
    # compensate with more substeps. Too few substeps here would blow up
    # to NaN (a stability failure, not a gradient bug -- see integrate.py).
    model = RRIModel(grid, n_substeps_slope=20, n_substeps_river=20)

    # Use a shallow soil layer + a heavier storm than the smoke test above,
    # so storage actually rises past `da` somewhere in the domain and the
    # Manning surface-flow branch (the only place `ns_slope` enters the
    # RHS) is exercised -- otherwise its gradient is *correctly* zero.
    soildepth = torch.full_like(params.slope.soildepth, 0.05)
    gammaa = params.slope.gammaa
    da = soildepth * gammaa
    dm = 0.3 * da
    params.slope.soildepth = soildepth
    params.slope.da = da
    params.slope.dm = dm

    ns_slope = params.slope.ns.clone().requires_grad_(True)
    ka = params.slope.ka.clone().requires_grad_(True)
    ns_river = params.river.ns_river.clone().requires_grad_(True)
    ksv = params.ksv.clone().requires_grad_(True)

    params.slope.ns = ns_slope
    params.slope.ka = ka
    params.river.ns_river = ns_river
    params.ksv = ksv

    T = 30
    rain = make_rain_pulse(T, grid.ny, grid.nx, peak_mm_per_hr=80.0, start=2, duration=10)
    dt = 300.0
    out = model.simulate(rain, params, dt, outlet_riv_index=outlet)
    loss = (out["qr_outlet"] ** 2).sum()
    loss.backward()

    for name, p in [("ns_slope", ns_slope), ("ka", ka), ("ns_river", ns_river), ("ksv", ksv)]:
        assert p.grad is not None, f"no gradient reached {name}"
        assert torch.isfinite(p.grad).all(), f"non-finite gradient for {name}"
        gnorm = p.grad.norm().item()
        print(f"[grad] |d(loss)/d({name})| = {gnorm:.6e}")
        assert gnorm > 0.0, f"zero gradient for {name}"


if __name__ == "__main__":
    grid, params, model, rain, dt, outlet, out = test_forward_runs()
    test_mass_balance(grid, params, model, rain, dt, outlet)
    test_gradients()
    print("\nAll checks passed.")
