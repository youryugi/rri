# rri_torch — a differentiable core of the RRI model

This is a from-scratch PyTorch reimplementation of the *core hydrology* of
the RRI (Rainfall-Runoff-Inundation) model whose Fortran/iRIC source lives
in `../src/RRI`. It is built for **gradient-based parameter calibration**:
every state update is a plain tensor operation, so you can backpropagate a
loss (e.g. MSE against an observed hydrograph) straight through the whole
time-stepping simulation to the physical parameters (Manning's roughness,
hydraulic conductivity, Green-Ampt infiltration parameters, ...).

It is **not** a drop-in replacement for the Fortran solver. It reproduces
the governing equations for the pieces listed below, using a different
(fixed-step) time integrator, and it does not implement dams, sediment
transport, levee breaks, diversions, groundwater exchange, or boundary
condition files. See "Scope and deviations" below before comparing results
against a real RRI run.

There is no bundled example DEM/rainfall/discharge dataset (the project
had none available). `examples/synthetic_basin.py` builds a small
synthetic V-shaped catchment instead — a straight river channel down the
middle draining two symmetric hillslopes — used by the tests and the demos.

## What's implemented, and where it comes from

| Module | Fortran reference | What it does |
|---|---|---|
| `geometry.py` | `RRI_Sub.f90` (`sub_slo_ij2idx`, `sub_riv_ij2idx`, D8 search) | Builds the domain/river masks and the river network's downstream tree from a D8 direction grid. Non-differentiable, built once. |
| `hydraulics.py` | `RRI_Slope.f90` (`hq`, `h2lev`), `RRI_Riv.f90` (`hq_riv`) | The tri-linear two-layer (unsaturated / saturated subsurface / Manning surface) storage-discharge law for slope cells, and Manning's equation for a rectangular river channel. |
| `slope.py` | `RRI_Slope.f90` (`funcs`, `qs_calc_do`) | The 8-direction diffusive-wave slope routing, vectorized as 4 half-edge flux tensors (E/S/SE/SW) plus their mirrored inflow, instead of the Fortran's explicit neighbour-index arrays. |
| `river.py` | `RRI_Riv.f90` (`funcr`, `qr_calc`) | Diffusive-wave routing along the river's D8 tree, using `index_add` scatter-add for conservative flux accumulation. |
| `exchange.py` | `RRI_RivSlo.f90` (`funcrs_dt`) | The 4-case weir-type river<->slope exchange. |
| `processes.py` | `RRI_Infilt.f90`, `RRI_Evp.f90` | Green-Ampt infiltration and simple potential-ET extraction. |
| `integrate.py` | `RRI_Mod2.f90`/`RRI.f90` (`runge_mod`, the adaptive stepping in the main loop) | Fixed-step RK4 integration (the default, for gradient-based work -- see below), plus a real adaptive Cash-Karp RKF45 integrator matching RRI's own step-size control exactly, for validation runs. |
| `model.py` | `RRI.f90` main loop | Assembles the above with the same operator-splitting order as the reference: river routing -> slope routing -> ET -> river/slope exchange -> infiltration -> drain outlet cells. |

`reference/` is a second, independent (no `rri_torch`/PyTorch import)
transcription of the same subroutines as plain Python/NumPy loops rather
than vectorized tensor ops, plus the real adaptive Cash-Karp RKF45
stepping -- see "Validation against an independent reference port" below.

## Scope and deviations from the reference Fortran

These are deliberate choices to keep a first differentiable version
tractable; each is a plausible place to extend later.

- **Not implemented at all:** groundwater exchange (`RRI_GW.f90`), dams
  (`RRI_Dam.f90`), diversions (`RRI_Div.f90`), levee breaks
  (`RRI_Break.f90`), sediment transport (`RRI_Sediment.f90`,
  `RRI_Sed2.f90`), non-rectangular channel cross-sections (`sec_map`
  lookup tables in `RRI_Section.f90`), the single-direction kinematic-wave
  slope mode (`eight_dir=0`/`dif=0`), and all boundary-condition files
  (discharge/water-level time series). Every slope cell is treated with
  the 8-direction diffusive scheme; every river reach is a plain
  rectangular channel.
- **Fixed-step RK4 by default, real adaptive RKF45 available.** The
  reference integrates both slope and river ODEs with an embedded,
  error-controlled Runge-Kutta-Fehlberg 4(5) scheme that shrinks its step
  until a local error estimate is satisfied. An adaptive step count that
  depends on the (learnable) parameters is awkward for a gradient-based
  workflow: it changes the computational graph's shape from run to run
  and can put kinks in the loss surface at step-accept/reject boundaries.
  `integrate.integrate_fixed` (the default -- `n_substeps_slope` /
  `n_substeps_river`) sidesteps this entirely; you're responsible for
  picking substep counts large enough (as you would size `dt`/`dt_riv` in
  the original model). **This is a real numerical stability tradeoff, not
  just an accuracy one:** with too few substeps the explicit scheme can
  diverge to NaN outright, especially with thin soil layers / large
  storms, or a real-world basin with a much stiffer reach than a
  synthetic test catchment ever exercises (see
  `tests/test_forward_and_grad.py::test_gradients` for a small-scale
  example needing 20 substeps where 4 blows up, and
  `differentiable/HANDOFF.md` section 6a for a real 15,751 km^2 basin
  needing 200 substeps on its stiffest reach). If you see NaNs, increase
  `n_substeps_*` before suspecting anything else. For validation/forward-
  only work where matching the reference's real step-size control matters
  more than a smooth loss surface, `integrate.integrate_adaptive` (wire
  up via `RRIModel(..., adaptive=True)`) implements RRI's actual Cash-Karp
  scheme -- gradients still flow through it (see its docstring for the
  autograd caveat), but prefer the fixed-step integrator for calibration.
- **The "don't overshoot" correction in the river/slope exchange** is an
  Fortran `do ... exit` loop of up to 10 iterations per cell; we replace
  it with a fixed `N_CORRECTION = 6` unrolled passes (a no-op past
  convergence), since a data-dependent iteration count has the same
  autograd-unfriendliness as the adaptive time-stepping above.
- **The overtopping-to-neighbours spreading** (`river_overtop_neighbor_switch`)
  and the channel-capacity/blockage extensions in `qr_calc` are optional,
  off-by-default experimental features in the reference; we don't
  implement them at all.
- **The Green-Ampt minimum-wetness cutoff is applied symmetrically.** In
  `RRI_Slope.f90::hq`, the `dh > 0 .and. h <= min_hs` cutoff that
  suppresses lateral subsurface flow below a wetness threshold is only
  ever true when the *current loop cell* is the uphill/outflow side of an
  edge; the call made for the downhill/inflow side is passed a negative
  `dh`, so the same check can never trigger there, even when the neighbor
  supplying the water is the one that's too dry. We apply the cutoff using
  whichever cell is actually supplying the water on each edge (see
  `hydraulics.hq_slope` and its call sites in `slope.py`), which we believe
  is the physically intended behaviour.
- **The outlet boundary condition is a free-flow assumption.** A river
  cell with no downstream cell drains against a synthetic bed slope
  (`RiverParams.outlet_slope`, default 1e-3) rather than an explicit
  water-level/discharge boundary file.
- **Safe square roots.** Every `sqrt` of a flow-driving depth/gradient is
  floored at a small epsilon (`tensor_ops.safe_sqrt`) rather than exactly
  zero. `sqrt(x)` has an infinite derivative at `x=0`; even where a branch
  is masked out by `torch.where` downstream, autograd still evaluates that
  branch's local gradient and multiplies it by a zero mask, and `0 * inf`
  is `NaN`. This is a purely numerical fix (negligible forward bias) but
  worth knowing about if you add new flux formulas.

## Validation against an independent reference port

There is no Intel Fortran compiler, `iriclib`, iRIC installation, or
sample `.cgn` project file available in this environment (see "Validating
against the Fortran model" below), so the compiled `install/rri.exe`
cannot actually be run here. As the next best thing, `reference/` is a
second, completely independent implementation of the same equations:
plain Python + NumPy, no `rri_torch`/PyTorch import, explicit per-cell
loops and if/elif branches instead of vectorized tensor ops (so it reads
as a much more literal transcription of the Fortran), and it also
implements the *real* adaptive Cash-Karp RKF45 stepping (the exact
coefficients from `RRI_Mod2.f90`), not just a fixed-step stand-in.

`tests/test_reference_agreement.py` checks two things against it, on the
same synthetic catchment:

1. **At matching (dt, n_substeps), `rri_torch` and `reference` agree to
   ~1e-17 relative error.** Since one is vectorized tensor ops and the
   other is scalar Python loops, this is strong evidence the vectorized
   translation is arithmetically faithful to the equations, not just
   "close."
2. **`rri_torch`'s fixed-step result converges to `reference`'s adaptive
   RKF45 result as `n_substeps` increases** (0.56% -> 0.017% -> 0.0008%
   -> 0.0002% at n_substeps = 4, 8, 16, 32), confirming the fixed-step
   scheme is a valid discretization of the same adaptive-RK physics RRI
   actually uses, not a different approximation that happens to look
   reasonable.
3. **A small explicit river confluence (two tributaries merging into one
   cell) agrees exactly** between `river_rhs`'s vectorized `index_add`
   inflow accumulation and `reference`'s plain dict `+=` loop
   (`test_river_confluence`) -- the synthetic catchment above is a single
   straight channel with no confluence, so this exercises a code path
   nothing else here does.

Separately, `differentiable/HANDOFF.md` (section 7) documents a real-data
validation against the actual compiled Fortran binary on a real 15,751
km^2 basin (not just the synthetic catchment), which is how two real
formula bugs were actually found and fixed -- worth reading if you're
about to trust this model's output against real observations.

**A finding worth knowing if you reuse `reference`'s adaptive mode on your
own catchment:** RRI's `eps` adaptive-step tolerance is an *absolute*
error bound in storage units (metres of depth for slope, and depth/volume
for river), not relative to the state's own magnitude. The Fortran
default (`eps = 0.01`, i.e. 1 cm) is sized for realistic river depths of
decimetres to metres. Our synthetic catchment's river only ever reaches a
few centimetres deep, so at the default `eps` the adaptive integrator
accepted just 1-4 giant substeps per 300 s window and disagreed with the
(already fully self-converged) fixed-step result by ~20% -- looking, at
first, exactly like a translation bug. Tightening `eps` from 0.01 to 1e-7
made that gap shrink smoothly to ~1e-6 with no floor, which is what
confirmed it was a tolerance/scenario mismatch rather than an error in
either implementation. If you point `reference`'s adaptive mode at a
catchment with realistically-sized depths, RRI's own `eps` convention
should be fine; for anything much smaller (as any minimal synthetic test
tends to be), tighten it accordingly.

## Quick start

```bash
cd differentiable
python -m venv .venv
./.venv/Scripts/python.exe -m pip install numpy torch --index-url https://download.pytorch.org/whl/cpu
./.venv/Scripts/python.exe tests/test_forward_and_grad.py         # forward + mass-balance + gradient sanity checks
./.venv/Scripts/python.exe tests/test_reference_agreement.py      # cross-check against the independent reference port
./.venv/Scripts/python.exe examples/calibrate_example.py          # gradient-based calibration demo (a few minutes on CPU)
```

## Using it against your own catchment

1. Build a `Grid` with `Grid.build(...)`: a `domain` mask, ground elevation
   `zb`, a D8 `dir` grid (for tracing the river network only), and a
   `riv_mask`. This mirrors what RRI's `.dir`/`.riv` input rasters encode.
2. Build `SlopeParams` / `RiverParams` / `RRIParams` — these are the
   learnable physical parameters (Manning's n, conductivities, Green-Ampt
   parameters, ...). Wrap whichever ones you want to calibrate as
   `torch.nn.Parameter` (or a reparametrized `exp(log_param)` as in
   `examples/calibrate_example.py` to keep them positive).
3. Run `RRIModel(grid).simulate(rain_seq, params, dt, outlet_riv_index=...)`
   and backpropagate a loss computed from `out["qr_outlet"]` (or any other
   state you record with `record_full_state=True`).

`rain_seq` is a `(T, ny, nx)` tensor of rainfall rate in m/s — RRI's
rainfall input files are typically mm/h; convert with `mm_per_hr * 1e-3 / 3600`.

## Validating against the Fortran model

This is "RRI on iRIC": `main.f90` requires a `.cgn` (CGNS) project file
that is normally authored and exported by the iRIC GUI, and
`install/rri.exe` (a working precompiled binary *is* present in this
repository) reads it through iRIC's proprietary `cg_iric_open` API. There
is no iRIC installation, no sample `.cgn` file, and no `iriclib` source in
this environment, so running the actual compiled solver on a matching
test case isn't possible here -- see "Validation against an independent
reference port" above for what we did instead (an independent loop-based
Python transcription of the same subroutines, including the real adaptive
RKF45 stepping, agreeing with `rri_torch` to ~1e-17 at matching step
sizes and converging to it as substeps increase).

If you do have a `.dir`/`.riv`/rainfall/discharge dataset the Fortran
solver already runs on (or can get one from iRIC), the most useful next
validation would be:

1. Feed it through `Grid.build` and the same parameters, run
   `RRIModel.simulate` with a fixed, small `dt` (well below the Fortran's
   typical `ddt_min`), and compare the outlet hydrograph against the real
   `rri.exe` run.
2. Expect small differences from: fixed vs. adaptive time-stepping, the
   symmetric min-wetness cutoff fix noted above, and the fixed-iteration
   exchange correction. If those differences are larger than expected,
   that's a real bug to chase down, not just "the two are different
   models."
