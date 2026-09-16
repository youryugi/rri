"""Cross-checks the vectorized, differentiable `rri_torch` model against
`reference`: an independent, unvectorized, loop-based Python port of the
same Fortran subroutines that has zero import dependency on `rri_torch`
or PyTorch (see `reference/*_ref.py`).

Two separate claims are checked here:

1. **The vectorized torch translation is arithmetically correct.**
   Run both models with the identical fixed-step RK4 scheme at the same
   dt/substep count -- they should agree to near machine precision, since
   at that point they're evaluating the exact same equations, just with
   torch tensor ops vs. plain Python loops.

2. **A fixed step is a reasonable stand-in for RRI's real adaptive
   RKF45 stepping.** The reference also implements the actual embedded
   Runge-Kutta-Fehlberg 4(5) (Cash-Karp) adaptive scheme RRI.f90 uses
   (see reference/integrate_ref.py), which is the closest thing to "what
   the compiled Fortran solver would produce" available without the
   actual binary + a CGNS project file (there is no sample dataset and no
   iRIC installation in this environment -- see the README). We show the
   torch fixed-step result converges toward that adaptive reference as
   `n_substeps` increases.

Run with:  .venv/Scripts/python.exe tests/test_reference_agreement.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch

from rri_torch.model import RRIModel
from reference.scenario_ref import build_synthetic_basin_ref
from reference.model_ref import RefModel
from examples.synthetic_basin import build_synthetic_basin, make_rain_pulse, DTYPE


def _rain_numpy(rain_t: torch.Tensor) -> np.ndarray:
    return rain_t.numpy()


def test_fixed_step_matches_reference(n_substeps=8, T=20, dt=300.0):
    grid, params, outlet = build_synthetic_basin()
    grid_ref, params_ref, outlet_key = build_synthetic_basin_ref()

    model = RRIModel(grid, n_substeps_slope=n_substeps, n_substeps_river=n_substeps)
    ref_model = RefModel(grid_ref, mode="fixed", n_substeps_slope=n_substeps, n_substeps_river=n_substeps)

    rain = make_rain_pulse(T, grid.ny, grid.nx)
    rain_np = _rain_numpy(rain)

    with torch.no_grad():
        out = model.simulate(rain, params, dt, outlet_riv_index=outlet, record_full_state=True)
    out_ref = ref_model.simulate(rain_np, params_ref, dt, outlet_key=outlet_key)

    qr_torch = out["qr_outlet"].numpy()
    qr_ref = out_ref["qr_outlet"]
    max_abs_diff = np.max(np.abs(qr_torch - qr_ref))
    scale = max(np.max(np.abs(qr_ref)), 1e-12)
    rel_diff = max_abs_diff / scale
    print(f"[fixed-step agreement] max|qr_torch - qr_ref| = {max_abs_diff:.3e}  "
          f"(relative to peak {scale:.3e}: {rel_diff:.3e})")
    assert rel_diff < 1e-6, "vectorized torch model disagrees with the loop-based reference at identical (dt, n_substeps)"

    hs_torch = out["hs"].numpy()
    hs_ref = out_ref["hs"]
    hs_diff = np.max(np.abs(hs_torch - hs_ref))
    print(f"[fixed-step agreement] max|hs_torch - hs_ref| over full trajectory = {hs_diff:.3e}")
    assert hs_diff < 1e-8, "slope-state trajectories disagree between torch and reference"


def test_converges_to_adaptive_reference(T=20, dt=300.0):
    """NOTE on `eps`: RRI's adaptive step-size control accepts a step
    when the estimated per-cell error is below `eps` *in absolute storage
    units* (metres of depth / m^3 of volume), not relative to the state's
    own magnitude -- see RRI.f90's `errmax = maxval(hs_err)/eps` and the
    equivalent for the river. The Fortran default (eps=0.01, i.e. 1 cm) is
    sized for realistic river depths of decimetres to metres. Our
    synthetic catchment's river only ever reaches a few *centimetres*
    deep, so that default tolerance is far too loose here -- accepting
    1-4 giant substeps per 300 s window instead of resolving anything.
    That was verified directly: shrinking `eps` from 0.01 down to 1e-7
    drives the adaptive result smoothly from a ~20% mismatch against the
    (already fully self-converged) fixed-step result down to ~1e-6, with
    no floor -- i.e. both integrators are solving the same equations
    correctly, and the gap at eps=0.01 was a tolerance/scenario mismatch,
    not a bug. We use eps=1e-6 below, appropriate for this catchment's
    scale; a real catchment with realistic depths should use RRI's own
    eps convention instead.
    """
    grid, params, outlet = build_synthetic_basin()
    grid_ref, params_ref, outlet_key = build_synthetic_basin_ref()

    rain = make_rain_pulse(T, grid.ny, grid.nx)
    rain_np = _rain_numpy(rain)

    ref_adaptive = RefModel(grid_ref, mode="adaptive", eps=1e-6, ddt_min_riv=1e-4, ddt_min_slo=1e-4)
    out_adaptive = ref_adaptive.simulate(rain_np, params_ref, dt, outlet_key=outlet_key)
    qr_adaptive = out_adaptive["qr_outlet"]
    scale = max(np.max(np.abs(qr_adaptive)), 1e-12)

    print("[convergence to adaptive RKF45 reference, eps=1e-6]")
    errors = []
    # n_substeps=2 is excluded: it's coarse enough to diverge to NaN for
    # this scenario, same fixed-step-stability caveat as noted in
    # rri_torch/integrate.py -- not a reference-agreement issue.
    for n_substeps in (4, 8, 16, 32):
        model = RRIModel(grid, n_substeps_slope=n_substeps, n_substeps_river=n_substeps)
        with torch.no_grad():
            out = model.simulate(rain, params, dt, outlet_riv_index=outlet)
        qr_torch = out["qr_outlet"].numpy()
        rel_err = np.max(np.abs(qr_torch - qr_adaptive)) / scale
        errors.append(rel_err)
        print(f"  n_substeps={n_substeps:3d}  rel_err_vs_adaptive={rel_err:.4%}")

    assert errors[-1] < errors[0], "fixed-step result should get closer to the adaptive reference as substeps increase"
    assert errors[-1] < 0.01, f"even at n_substeps=32 the gap to the (tight-tolerance) adaptive reference is large: {errors[-1]:.2%}"


if __name__ == "__main__":
    test_fixed_step_matches_reference()
    test_converges_to_adaptive_reference()
    print("\nAll reference-agreement checks passed.")
