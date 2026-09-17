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
    T = 80
    rain = make_rain_pulse(T, grid.ny, grid.nx, start=4, duration=12)
    dt = 150.0  # 2.5 minutes -- see test_mass_balance for why not 300s

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
    #
    # Why dt=150s and not the old 300s: the river<->slope exchange, ET, and
    # sink-drain steps run once per *outer* dt regardless of substep count
    # (model.py's operator-splitting order matches RRI.f90's own), so their
    # contribution to the splitting error is O(dt), not shrinkable by
    # raising n_substeps_slope/river (verified directly: ns 4->16 at a
    # fixed dt=300s left the error flat at ~8.4%, while halving dt at a
    # fixed ns=4 dropped it to <1%). This got a lot more visible after
    # fixing exchange.py's missing "both banks" factor of 2 (see
    # HANDOFF.md section 6a) -- exchange moving ~2x more water per outer
    # step doubles this splitting error's absolute size at the same dt.
    # 150s keeps this smoke test's error comfortably under the 5% bar
    # without having to also raise n_substeps (which wouldn't help).
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


def _build_stiff_gradient_scenario():
    """Shallow soil layer + heavier storm than the smoke test above, so
    storage actually rises past `da` somewhere in the domain and the
    Manning surface-flow branch (the only place `ns_slope` enters the
    RHS) is exercised -- otherwise its gradient is *correctly* zero. This
    is also numerically stiffer than the smoke test, which is the point:
    it exercises the same regime that needs either enough fixed substeps
    or (see `test_gradients_adaptive`) real adaptive stepping to avoid
    diverging to NaN.
    """
    grid, params, outlet = build_synthetic_basin()
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

    grad_params = {"ns_slope": ns_slope, "ka": ka, "ns_river": ns_river, "ksv": ksv}
    return grid, params, outlet, grad_params


def _check_gradients(model, grid, params, outlet, grad_params, tag, use_checkpointing=False):
    T = 30
    rain = make_rain_pulse(T, grid.ny, grid.nx, peak_mm_per_hr=80.0, start=2, duration=10)
    dt = 300.0
    out = model.simulate(rain, params, dt, outlet_riv_index=outlet, use_checkpointing=use_checkpointing)
    loss = (out["qr_outlet"] ** 2).sum()
    loss.backward()

    for name, p in grad_params.items():
        assert p.grad is not None, f"no gradient reached {name}"
        assert torch.isfinite(p.grad).all(), f"non-finite gradient for {name}"
        gnorm = p.grad.norm().item()
        print(f"[grad{tag}] |d(loss)/d({name})| = {gnorm:.6e}")
        assert gnorm > 0.0, f"zero gradient for {name}"


def test_gradients():
    # This scenario is numerically stiffer than the smoke test above -- the
    # reference Fortran handles that by shrinking its adaptive RK step; we
    # use a fixed step, so we compensate with more substeps. Too few
    # substeps here would blow up to NaN (a stability failure, not a
    # gradient bug -- see integrate.py).
    grid, params, outlet, grad_params = _build_stiff_gradient_scenario()
    model = RRIModel(grid, n_substeps_slope=20, n_substeps_river=20)
    _check_gradients(model, grid, params, outlet, grad_params, tag="")


def test_gradients_adaptive():
    """Same stiff scenario, but through `integrate_adaptive` (see
    HANDOFF.md section 6a and integrate.py's module docstring): confirms
    gradients still flow correctly through the real Cash-Karp RKF45 path,
    not just fixed-step RK4. `eps`/`ddt_min` are tightened from RRI's own
    defaults for the same reason `test_reference_agreement.py` tightens
    them -- this synthetic catchment's depths are centimetre-scale, far
    below what the 1cm-absolute default tolerance was sized for.
    """
    grid, params, outlet, grad_params = _build_stiff_gradient_scenario()
    model = RRIModel(grid, adaptive=True, eps=1e-6, ddt_min_slope=1e-4, ddt_min_river=1e-4)
    _check_gradients(model, grid, params, outlet, grad_params, tag="-adaptive")


def test_gradients_checkpointed_match_uncheckpointed():
    """`simulate(..., use_checkpointing=True)` trades ~2x forward compute
    (each step's internal substeps are recomputed once during backward)
    for O(n_substeps)->O(1) memory per step -- confirmed necessary on a
    real (not toy) grid: backpropagating through just 200 un-checkpointed
    outer steps at n_substeps_slope=40 on the real Solo scenario's
    18582-cell domain (see HANDOFF.md's calibration-plan notes) OOMs a
    32GB GPU outright. Checkpointing is mathematically exact (not an
    approximation), so this checks the gradients it produces are
    *identical* to the un-checkpointed path on this small synthetic
    scenario, not just "close enough".
    """
    grid, params, outlet, grad_params = _build_stiff_gradient_scenario()
    model = RRIModel(grid, n_substeps_slope=20, n_substeps_river=20)
    T = 30
    rain = make_rain_pulse(T, grid.ny, grid.nx, peak_mm_per_hr=80.0, start=2, duration=10)
    dt = 300.0
    out = model.simulate(rain, params, dt, outlet_riv_index=outlet, use_checkpointing=True)
    loss = (out["qr_outlet"] ** 2).sum()
    loss.backward()

    # exact values from test_gradients() (same scenario, same seed-free
    # deterministic setup) -- checkpointing must reproduce them exactly.
    expected = {"ns_slope": 1.172758e+00, "ka": 2.216010e-01, "ns_river": 1.956847e+02, "ksv": 4.376833e+05}
    for name, p in grad_params.items():
        gnorm = p.grad.norm().item()
        print(f"[grad-checkpointed] |d(loss)/d({name})| = {gnorm:.6e}")
        rel_err = abs(gnorm - expected[name]) / expected[name]
        assert rel_err < 1e-6, f"checkpointed gradient for {name} disagrees with un-checkpointed: {gnorm} vs {expected[name]}"


def test_checkpoint_substeps_match_uncheckpointed():
    """`RRIModel(checkpoint_substeps=True)` nests a second level of
    checkpointing inside `integrate_fixed`'s RK4 substep loop (see its
    docstring) -- needed because on a real grid, `n_substeps_{slope,river}`
    is large enough that even *one* outer step's local backward graph
    (all substeps unrolled) is too big to hold, independent of outer-step
    count (confirmed OOMing at 25+ GB on the real Solo scenario at just
    T=5 outer steps -- see HANDOFF.md). Like outer-step checkpointing,
    this is mathematically exact, not an approximation: checked here
    against the same hardcoded expected gradients as
    `test_gradients_checkpointed_match_uncheckpointed`, with *both*
    checkpoint levels turned on together (the real usage pattern).
    """
    grid, params, outlet, grad_params = _build_stiff_gradient_scenario()
    model = RRIModel(grid, n_substeps_slope=20, n_substeps_river=20, checkpoint_substeps=True)
    _check_gradients(model, grid, params, outlet, grad_params, tag="-nested-checkpoint", use_checkpointing=True)

    expected = {"ns_slope": 1.172758e+00, "ka": 2.216010e-01, "ns_river": 1.956847e+02, "ksv": 4.376833e+05}
    for name, p in grad_params.items():
        gnorm = p.grad.norm().item()
        rel_err = abs(gnorm - expected[name]) / expected[name]
        assert rel_err < 1e-6, f"nested-checkpoint gradient for {name} disagrees: {gnorm} vs {expected[name]}"


def test_checkpoint_substeps_rejects_track_qr_avg():
    grid, params, outlet = build_synthetic_basin()
    try:
        RRIModel(grid, checkpoint_substeps=True, track_qr_avg=True)
        assert False, "expected ValueError for checkpoint_substeps + track_qr_avg"
    except ValueError:
        pass


if __name__ == "__main__":
    grid, params, model, rain, dt, outlet, out = test_forward_runs()
    test_mass_balance(grid, params, model, rain, dt, outlet)
    test_gradients()
    test_gradients_adaptive()
    test_gradients_checkpointed_match_uncheckpointed()
    test_checkpoint_substeps_match_uncheckpointed()
    test_checkpoint_substeps_rejects_track_qr_avg()
    print("\nAll checks passed.")
