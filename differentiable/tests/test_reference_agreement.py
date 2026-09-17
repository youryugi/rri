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
from rri_torch.geometry import Grid
from rri_torch.river import RiverParams, river_rhs
from reference.scenario_ref import build_synthetic_basin_ref
from reference.model_ref import RefModel
from reference.geometry_ref import build_ref_grid
from reference.physics_ref import river_rhs_ref
from reference.hydraulics_ref import hq_river as hq_river_ref
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


def test_adaptive_matches_reference(T=20, dt=300.0):
    """`rri_torch.integrate.integrate_adaptive` (the real Cash-Karp RKF45
    port, see HANDOFF.md section 6a) should agree with the independent
    loop-based reference's own adaptive integrator to near machine
    precision when run with the exact same (eps, ddt_min) -- unlike
    `test_converges_to_adaptive_reference` below, this isn't a
    convergence claim, it's the same arithmetic-correctness claim as
    `test_fixed_step_matches_reference` but for the adaptive path.
    """
    grid, params, outlet = build_synthetic_basin()
    grid_ref, params_ref, outlet_key = build_synthetic_basin_ref()

    eps, ddt_min = 1e-6, 1e-4
    model = RRIModel(grid, adaptive=True, eps=eps, ddt_min_slope=ddt_min, ddt_min_river=ddt_min)
    ref_model = RefModel(grid_ref, mode="adaptive", eps=eps, ddt_min_riv=ddt_min, ddt_min_slo=ddt_min)

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
    print(f"[adaptive agreement] max|qr_torch - qr_ref| = {max_abs_diff:.3e}  "
          f"(relative to peak {scale:.3e}: {rel_diff:.3e})")
    assert rel_diff < 1e-6, "adaptive torch integrator disagrees with the loop-based adaptive reference"

    hs_torch = out["hs"].numpy()
    hs_ref = out_ref["hs"]
    hs_diff = np.max(np.abs(hs_torch - hs_ref))
    print(f"[adaptive agreement] max|hs_torch - hs_ref| over full trajectory = {hs_diff:.3e}")
    assert hs_diff < 1e-8, "slope-state trajectories disagree between torch and reference under adaptive stepping"


def test_river_confluence():
    """`river_rhs`'s inflow accumulation (`inflow.index_add(0, down[has_down],
    qr[has_down])`, river.py) has never actually been exercised at a real
    confluence by any existing test: `examples.synthetic_basin` is a single
    straight channel with no two reaches merging into one cell, so a
    duplicate-index bug in the vectorized `index_add` scatter could have
    gone undetected. This builds a small explicit Y-shaped network --
    two tributaries (0,0) and (0,2) both draining into a confluence cell
    (1,1), which continues down to an outlet -- and checks:

    1. `river_rhs` (torch, vectorized `index_add`) agrees with
       `river_rhs_ref` (plain Python dict accumulation, `inflow[down] +=
       qr[(i,j)]`) to near machine precision.
    2. The confluence cell's inflow independently hand-computed from the
       two tributaries' own `hq_river` outputs (via the reference's own
       hydraulics function, so this isn't circular) matches what both
       implementations report.
    """
    ny, nx = 4, 3
    domain = np.ones((ny, nx), dtype=bool)
    riv_mask = np.zeros((ny, nx), dtype=bool)
    riv_mask[0, 0] = riv_mask[0, 2] = True   # two tributary heads
    riv_mask[1, 1] = True                     # confluence
    riv_mask[2, 1] = riv_mask[3, 1] = True    # mainstem to the outlet

    dir_grid = np.zeros((ny, nx), dtype=np.int64)
    dir_grid[0, 0] = 2   # SE -> (1,1)
    dir_grid[0, 2] = 8   # SW -> (1,1)
    dir_grid[1, 1] = 4   # S  -> (2,1)
    dir_grid[2, 1] = 4   # S  -> (3,1)
    dir_grid[3, 1] = 0   # outlet

    dx = dy = 100.0
    zb = 10.0 - 2.0 * np.repeat(np.arange(ny)[:, None], nx, axis=1).astype(np.float64)
    len_riv = np.where(riv_mask, dx, 0.0)

    grid = Grid.build(domain=domain, zb=zb, dir_grid=dir_grid, riv_mask=riv_mask,
                       len_riv=len_riv, dx=dx, dy=dy, dtype=DTYPE)
    grid_ref = build_ref_grid(domain=domain, zb=zb, dir_grid=dir_grid, riv_mask=riv_mask,
                               len_riv=len_riv, dx=dx, dy=dy)

    width, depth, ns_river, height, outlet_slope = 5.0, 1.0, 0.03, 0.0, 1e-3
    n_riv = grid.n_riv
    river_params = RiverParams(
        width=torch.full((n_riv,), width, dtype=DTYPE),
        depth=torch.full((n_riv,), depth, dtype=DTYPE),
        ns_river=torch.tensor(ns_river, dtype=DTYPE),
        height=torch.full((n_riv,), height, dtype=DTYPE),
        outlet_slope=outlet_slope,
    )
    river_params_ref = dict(
        width={k: width for k in grid_ref.down_riv},
        depth={k: depth for k in grid_ref.down_riv},
        height={k: height for k in grid_ref.down_riv},
        ns_river=ns_river,
        outlet_slope=outlet_slope,
    )

    # distinct depths everywhere a river cell exists, so no cancellation
    # accidentally hides a confluence-summation bug.
    hr_by_ij = {(0, 0): 0.50, (0, 2): 0.30, (1, 1): 0.10, (2, 1): 0.04, (3, 1): 0.0}
    hr_torch = torch.zeros(n_riv, dtype=DTYPE)
    for (i, j), h in hr_by_ij.items():
        hr_torch[grid.riv_index_map[i, j]] = h
    hr_ref = dict(hr_by_ij)

    dhrdt_torch, qr_torch = river_rhs(hr_torch, grid, river_params)
    dhrdt_ref, qr_ref = river_rhs_ref(hr_ref, grid_ref, river_params_ref)

    def torch_at(field, i, j):
        return field[grid.riv_index_map[i, j]].item()

    max_diff = 0.0
    for (i, j) in hr_by_ij:
        d_t, d_r = torch_at(dhrdt_torch, i, j), dhrdt_ref[(i, j)]
        q_t, q_r = torch_at(qr_torch, i, j), qr_ref[(i, j)]
        max_diff = max(max_diff, abs(d_t - d_r), abs(q_t - q_r))
    print(f"[confluence] max|torch - ref| over dhr/dt and qr at all 5 cells = {max_diff:.3e}")
    assert max_diff < 1e-12, "vectorized river_rhs disagrees with the loop-based reference at a confluence"

    # Independent hand-check: the confluence cell (1,1) should receive
    # exactly qr(0,0)+qr(0,2) as inflow, computed via the reference's own
    # (already-verified-elsewhere) hq_river, not by re-deriving qr_ref above.
    zb_conf = zb[1, 1] - depth
    for (i, j) in [(0, 0), (0, 2)]:
        zb_p = zb[i, j] - depth
        dist = grid_ref.dis_riv[(i, j)]
        dh = ((zb_p + hr_by_ij[(i, j)]) - (zb_conf + hr_by_ij[(1, 1)])) / dist
        assert dh >= 0.0, "test setup should have flow running downhill into the confluence"
    q00 = hq_river_ref(hr_by_ij[(0, 0)], ((zb[0, 0] - depth + hr_by_ij[(0, 0)]) - (zb_conf + hr_by_ij[(1, 1)])) / grid_ref.dis_riv[(0, 0)], width, ns_river)
    q02 = hq_river_ref(hr_by_ij[(0, 2)], ((zb[0, 2] - depth + hr_by_ij[(0, 2)]) - (zb_conf + hr_by_ij[(1, 1)])) / grid_ref.dis_riv[(0, 2)], width, ns_river)
    channel_area_conf = width * grid_ref.len_riv_grid[1, 1]
    q_out_conf = torch_at(qr_torch, 1, 1)
    expected_dhrdt_conf = (-q_out_conf + q00 + q02) / channel_area_conf
    print(f"[confluence] hand-computed inflow = {q00 + q02:.6e}, "
          f"dhr/dt(confluence) expected={expected_dhrdt_conf:.6e} vs torch={torch_at(dhrdt_torch, 1, 1):.6e}")
    assert abs(torch_at(dhrdt_torch, 1, 1) - expected_dhrdt_conf) < 1e-12, \
        "confluence inflow does not equal the sum of both tributaries' independently hand-computed discharge"


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
    test_adaptive_matches_reference()
    test_river_confluence()
    test_converges_to_adaptive_reference()
    print("\nAll reference-agreement checks passed.")
