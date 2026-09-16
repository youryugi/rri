# Handoff: Differentiable RRI — status and next steps

Written for: continuing this work in a new session (e.g. on a GPU machine), so it
needs to stand alone without the prior conversation history.

## 1. What this project is

`differentiable/rri_torch/` is a from-scratch PyTorch reimplementation of the
*core hydrology* of the RRI (Rainfall-Runoff-Inundation) model — the Fortran/iRIC
source lives in `../src/RRI` of this repo. It's built so that a full rainfall-to-
discharge simulation is a plain sequence of differentiable tensor ops, enabling
gradient-based calibration of physical parameters (Manning's roughness, lateral
conductivity, Green-Ampt infiltration parameters, ...) by backpropagating a loss
against observed hydrographs — no black-box/gradient-free search needed.

Full design rationale, scope, and deliberate deviations from the Fortran are in
`differentiable/README.md`. Read that first if picking this up cold.

## 2. What's built and verified so far

- `rri_torch/`: slope 2D diffusive-wave routing, river 1D diffusive-wave routing
  (D8 tree), river↔slope weir exchange, Green-Ampt infiltration, simple ET,
  fixed-step RK4 time integration, all as vectorized PyTorch tensor ops.
- `reference/`: an independent, non-vectorized, loop-based Python/NumPy port of
  the same equations (no `rri_torch`/PyTorch import), including the *real*
  adaptive Cash-Karp RKF45 stepping (not just a fixed-step stand-in). Exists to
  cross-check `rri_torch`'s vectorized translation.
- `examples/synthetic_basin.py`: builds a small synthetic 12x7 V-shaped test
  catchment (no real data existed at that point in the project).
- `examples/calibrate_example.py`, `examples/compare_gradient_vs_gradientfree.py`:
  demonstrate gradient-based calibration and benchmark it against a gradient-free
  baseline (Nelder-Mead) on the synthetic catchment.
- Tests, all passing as of this handoff:
  - `tests/test_forward_and_grad.py`: forward run, mass balance (~3% residual,
    expected from operator splitting), and finite/non-zero gradients w.r.t.
    ns_slope, ka, ns_river, ksv.
  - `tests/test_reference_agreement.py`: `rri_torch` vs `reference` agree to
    ~1e-17 relative error at matching (dt, n_substeps); `rri_torch`'s fixed-step
    result converges to `reference`'s adaptive-RKF45 result as n_substeps
    increases (0.56% → 0.0002% at n_substeps = 4 → 32). **Important finding**:
    RRI's adaptive-step `eps` is an *absolute* tolerance (metres of storage),
    not relative — the Fortran default (`eps=0.01`) is sized for realistic river
    depths and is far too loose for a small synthetic catchment's centimetre-
    scale depths. Use a tighter `eps` (we used 1e-6) when comparing against a
    small/synthetic scenario. See the README for the full writeup.

Setup: `cd differentiable && python -m venv .venv && ./.venv/Scripts/python.exe -m pip install numpy torch --index-url https://download.pytorch.org/whl/cpu` (swap the index URL for a CUDA build on the GPU machine).

## 3. What we were doing when this session ended: real-data validation

Everything above was only validated against a synthetic catchment and against our
own independent reference port — never against real observed data, and never
against the actual compiled Fortran binary. We went looking for a real dataset to
close that gap.

### 3.1 Why the *iRIC* build of RRI (this repo's main src) couldn't be used

`src/RRI` + `main.f90` in this repo is "RRI on iRIC": it requires a `.cgn` (CGNS)
project file normally authored via the iRIC GUI, read through iRIC's proprietary
`cg_iric_open` API. No iRIC installation, no sample `.cgn` file, and no `iriclib`
source were available, so the precompiled `install/rri.exe` in this repo could not
actually be run here.

### 3.2 What we found instead: the official RRI-GUI distribution

The user downloaded the official RRI distribution from PWRI/ICHARM (requires
accepting a terms-of-use page with email/country at
`https://www.pwri.go.jp/icharm/research/rri/rri_contract_e.html` — already done)
to:

```
C:\Users\yang\Desktop\github\rri\RRI_1_4_2_7_GUI_Beta\
```

Inside it, two things turned out to be gold:

1. **`RRI-CUI/`** — the *classic, non-iRIC* RRI: plain-text ASCII I/O, no CGNS
   dependency at all. `RRI-CUI/bin/0_rri_1_4_2_6.exe` is a real, runnable
   Fortran binary.
2. **`RRI-CUI/Project/solo30s/`** — a **complete, ready-to-run real project**:
   the Solo River basin, Java, Indonesia (~15,751 km², 336×204 grid @ 30 arcsec,
   18,582 active cells). This is the basin used in RRI's own official tutorial
   (Sayama & Iwami, 2016, `RRI_CHA.pdf`, fetched during this session from
   `https://hywr.kuciv.kyoto-u.ac.jp/ihp/cha/tools/RRI/RRI_CHA.pdf`).

We **successfully ran the real compiled binary** in this environment:

```bash
cd "RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s"
cp "../../../../install/libiomp5md.dll" .   # fixes "error while loading shared libraries: libiomp5md.dll"
./0_rri_1_4_2_6.exe RRI_Input.txt            # exit code 0, ~360h/2160 steps in well under a minute
```

This produced real reference output in `solo30s/out/` — `hs_000001.out` ...
`hs_000096.out` (and same for `hr_`, `qr_`), plus `storage.dat`, one file per of
the 96 output steps (lasth=360h, outnum=96 → 3.75h between outputs).

**Also found real observed data** in `solo30s/obs/`: discharge and water-level
records at 6 gauge stations (`disc_jurug.data`, `disc_cepu.data`,
`disc_kajangan.data`, `disc_napel.data`, `disc_sekayu.data`, `disc_serenan.data`,
and matching `wl_*.data`), plus station grid coordinates in
`obs/location_solo_30s_all.txt` (`name row col`, e.g. `Jurug 117 79`). Format is
plain text, `iostat`-terminated, pairs of `<index> <value>` per line (see
`RRI-CUI/etc/evalHydro/evalHydro.f90` for the exact read pattern) — only 16
points per station (ties to rain.dat only having 16 daily blocks, see below), so
this is *not* a fine-grained continuous hydrograph; may be a peak/daily
comparison series rather than sub-daily. Not fully resolved — see open questions.

**Also found the original (non-iRIC) Fortran source**, cleaner than this repo's
heavily-modified iRIC version, at:

```
RRI_1_4_2_7_GUI_Beta/RRI-CUI/source/1.4.2.7/*.f90
```

Reading it clarified several details our `rri_torch` port either didn't handle or
handled differently from what a real lat-lon project needs:

| Detail | Real RRI behaviour (from `RRI.f90` v1.4.2.7) | Status in `rri_torch` |
|---|---|---|
| dx/dy for lat-lon (utm=0) grids | Computed **once for the whole domain** via the Hubeny geodesic-distance formula applied to the 4 bounding-box edges, then `dx = mean(N,S edge length)/nx`, `dy = mean(E,W edge length)/ny`. Not per-row/per-latitude. | **Not implemented** — `Grid.build` only takes a flat `dx`, `dy` in metres; caller must compute these themselves for lat-lon inputs. Need a small Hubeny-formula helper. |
| River reach length (no `len_riv` file given) | `len_riv = sqrt(dx * dy)` for **every** river cell, direction-independent. | Our synthetic basin passed `len_riv = dy` (a simplification for a north-south synthetic river). `Grid.build`'s own default is `dx`. **Should use `sqrt(dx*dy)` to match real RRI when no file is given.** |
| River width/depth (no file given) | `width = width_param_c * (acc * dx * dy * 1e-6) ** width_param_s`, same form for `depth` with its own params. Note the `acc * dx * dy * 1e-6` term is **drainage area in km²**, not raw accumulation cell count. | Not implemented as a formula — `RiverParams.width/depth` are plain per-cell tensors the caller must supply. Need to compute this from `acc` before building `RiverParams` for a real scenario. |
| River cell mask | `riv = 1` where `acc > riv_thresh` (strict `>`). | Matches what `Grid.build` expects the caller to pass as `riv_mask`; just need to compute it from `acc` the same way. |
| Rainfall file format (`rain.dat`) | Repeating blocks: header line `t_seconds  nx_rain  ny_rain`, then `ny_rain` rows of `nx_rain` values in mm/h. Converted to m/s via `/3600/1000`. Between blocks, rain is a **step function** (holds the last block's value) — confirmed from `RRI.f90`'s `t_rain(jtemp-1) < time+ddt <= t_rain(jtemp)` bin lookup, not interpolated. | Not implemented — need a small parser. Solo's `rain.dat` has exactly 16 blocks at t=0, 86400, 172800, ... (daily, in seconds), covering the full 360h = 15 days. |
| DEM/acc/dir ASCII format | Standard ESRI ASCII grid: `ncols`/`nrows`/`xllcorner`/`yllcorner`/`cellsize`/`NODATA_value` header, then rows of comma-space-separated values. | Not implemented — need a small parser (trivial). |
| This specific project's key parameters (`RRI_Input.txt`) | `utm=0`, `eight_dir=1` (diffusive, matches our implemented mode), `lasth=360`, `dt=600`, `dt_riv=60`, `ns_river=0.030`, `num_of_landuse=1`, `dif=1`, `ns_slope=0.400`, `soildepth=1.000`, `gammaa=0.475`, **`ksv=0.000` (kv), `faif=0.316` (Sf)**, **`ka=0.000`**, `gammam=0.000`, `beta=8.000`, `riv_thresh=100`, `width_param_c=5.000`, `width_param_s=0.350`, `depth_param_c=0.950`, `depth_param_s=0.200`, `height_param=0.000`, `height_limit_param=20.0`. Dam/diversion/boundary-files/land-use-file/PET/cross-section all **off (switch 0)**. | `ka=0` and `ksv=0` means this real case's subsurface/infiltration terms are numerically inert — it's effectively a pure surface-Manning + river test, which happens to line up exactly with `rri_torch`'s implemented scope (no GW, no dam, no boundary files, no land-use table, no cross-sections, no PET) — this real case is a very good fit for what's already built. |

The Solo project's `RRI_Input.txt` is at:
`RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s/RRI_Input.txt` (full contents were
read and are reproduced in the table above; re-read the file directly for exact
formatting/line numbers if writing a parser).

## 4. Open questions / where we got interrupted

The user cut off an `AskUserQuestion` about how to scope the real-data validation
run, given the real Solo basin is much bigger (18,582 active cells, 2160 time
steps) than anything `rri_torch` has been performance-tested on (the synthetic
catchment is ~60-80 cells). The options on the table were:

1. **Crop to a small sub-catchment** upstream of one gauge (e.g. Jurug, at grid
   row 117, col 79) and/or a short time window — fastest to get a first real
   comparison, but not a full-basin validation.
2. **Full basin, full 15-day event** — most convincing if it runs in reasonable
   time, but untested at this scale; may need performance work (this is
   presumably part of *why* the user is moving to a GPU machine).
3. **Full basin, but only the first 10-20 time steps** — cheap way to validate
   spatial correctness (compare `hs`/`hr` grids directly against the real
   Fortran's `out/hs_000001.out` etc.) without running the full event.

**This decision was never made — pick up here.** Given a GPU is now available,
option 2 (full basin, full event) is probably worth attempting directly, since
GPU throughput removes most of the original concern; fall back to option 1 or 3
if it's still too slow or something doesn't check out numerically.

## 5. Concrete next steps (in likely order)

1. Move/access `RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s/` (real input data
   + real Fortran output already generated in `out/`) alongside this repo on the
   new machine.
2. Write a small loader (new file, e.g. `differentiable/real_data/solo_river.py`):
   - ESRI ASCII grid parser (`ncols`/`nrows`/`xllcorner`/`yllcorner`/`cellsize`/
     `NODATA_value` header + comma-separated rows) for `topo/adem.txt`,
     `topo/acc_mod.txt`, `topo/dir_mod.txt`.
   - Hubeny-formula dx/dy computation from the bounding box (see table above;
     the exact subroutine is `hubeny_sub` in
     `RRI_1_4_2_7_GUI_Beta/RRI-CUI/source/1.4.2.7/RRI.f90`, read that for the
     precise formula before implementing).
   - Derive `riv_mask = acc > riv_thresh(=100)`, then `width`/`depth` via the
     power-law formula above, `len_riv = sqrt(dx*dy)` everywhere, `height = 0`.
   - `rain.dat` parser (repeating `t_seconds nx ny` + grid blocks; note Solo's
     rain grid extent might differ from the topo grid — check
     `xllcorner_rain`/`yllcorner_rain`/`cellsize_rain_x/y` in `RRI_Input.txt`
     against the topo grid's own corner/cellsize before assuming they're
     identical, and implement the nearest-cell index mapping RRI itself uses
     — see `rain_i`/`rain_j` computation in `RRI.f90` around the block that
     reads rain.dat).
   - `obs/disc_*.data` parser (plain `<index> <value>` pairs) and
     `obs/location_solo_30s_all.txt` parser (`name row col`).
3. Build a `Grid` + `RRIParams` from this real data (reuse
   `rri_torch.geometry.Grid.build`, `rri_torch.model.RRIParams`, etc. — the
   dataclasses are flexible enough, nothing in `rri_torch`'s core physics needs
   to change, only the scenario-construction script).
4. Decide the scope question from §4 and run `RRIModel.simulate` on it.
5. Compare against:
   - The real Fortran's `out/hs_*.out` / `hr_*.out` / `qr_*.out` grids directly
     (spatial-field agreement), and/or
   - The real observed discharge in `obs/disc_*.data` at the 6 gauge locations
     (the actual validation the user originally asked for).
6. Resolve the "16 points, what time spacing" question for the obs data if it
   matters for the comparison chosen (cross-reference against `rain.dat`'s 16
   daily blocks — they're very likely the same 16 daily timestamps, but confirm
   by reading `RRI-CUI/etc/calcHydro`/`evalHydro` source and/or the actual
   `disc_*.data` values against `hydro.txt`/`hydro_hr.txt` which already exist
   in the `solo30s` project folder as some kind of precomputed hydrograph
   output — worth reading those two files first, they may already contain a
   ready-made simulated-vs-observed comparison from ICHARM's own reference run).

## 6a. Update (2026-09-16 session, on the GPU machine): loader built, first
real-data run done, one open discrepancy found

Picked this up on the Linux GPU machine (RTX 5000 Ada, 32GB). The
`RRI_1_4_2_7_GUI_Beta` distribution (with `solo30s/` inside it) is now at
`/home/yang/github/rri/RRI_1_4_2_7_GUI_Beta/` in this repo (unzipped from
`RRI_1_4_2_7_GUI_Beta.zip`, checked into the working tree but not the
`differentiable/` subfolder itself).

**What's done:**

- `differentiable/real_data/solo_river.py`: the real-data loader described
  in section 5 is now written and working. ESRI ASCII grid parser, exact
  Hubeny geodesic dx/dy (GRS80 ellipsoid, `a=6378137`, `b=6356752.314` --
  transcribed directly from `hubeny_sub` in
  `RRI_1_4_2_7_GUI_Beta/RRI-CUI/source/1.4.2.7/RRI_Sub.f90` line ~499),
  `rain.dat` step-function parser with the exact `rain_i`/`rain_j` index
  mapping formula from `RRI.f90` lines 529-534, `width`/`depth` power-law
  formulas, `len_riv = sqrt(dx*dy)`, domain/sink masks (`zs > -100`,
  `dir in {0,-1}`), and obs/location file parsers. `build_solo_scenario()`
  returns a ready-to-use `(grid, params, outlet_riv_indices, rain)` tuple.
- **Validated the loader against the real compiled Fortran binary's own
  printed output** (see below): dx, dy, num_of_cell, and total area all
  matched to full double precision / exactly. This confirms the geometry
  half of the port (grid construction, Hubeny formula, domain/river
  masking) is bit-correct.
- Fixed a real GPU bug in `rri_torch/geometry.py`:
  `Grid.scatter_from_riv` created its output tensor with
  `torch.full(..., dtype=...)` but no `device=`, so it silently defaulted
  to CPU and crashed (`RuntimeError: ... same device ...`) the first time
  the model ran on CUDA. Fixed to `device=riv_field.device`. This bug is
  in the core library, not the new real-data code -- worth remembering if
  GPU runs on the synthetic scenario are ever tried too.
- `differentiable/examples/solo_validation.py`: runs `RRIModel` on the
  real Solo scenario for N steps, on CPU or CUDA, reports discharge at
  all 6 gauges, saves the series to a `.pt` file with `--save`.
- `differentiable/examples/compare_solo_hydro.py`: compares a saved run's
  Cepu discharge against `solo30s/hydro.txt` (see below -- this is RRI's
  own built-in point-hydrograph output, hourly, written directly by
  `RRI.f90` itself at the single station in `location.txt`, *not* derived
  from the `out/qr_*.out` spatial dumps -- search `hydro_file` in
  `RRI.f90` to see it's hardcoded). Computes NSE/RMSE/bias/peak-timing.

**Stability finding (expected, per section 6's gotcha, but the real
magnitude matters):** fixed-step RK4 with `n_substeps=20` (fine for the
synthetic catchment) diverges to NaN on the real Solo domain around
simulated hour 79 (right as the heaviest rain block, 8.33mm/h at t=48h,
is still being routed through the network -- the real terrain is Java
volcanic topography, much steeper/stiffer than the synthetic V-shaped
catchment). **`n_substeps=40` (both slope and river) ran the entire
360h/2160-step event with no NaN** in ~24 minutes on the RTX 5000 Ada.
Confirmed via a short bisection (20 fails ~step 475, 40 stable through
step 550 which spans the whole peak-intensity block) before committing
GPU time to the full run -- worth doing this bisection-on-a-short-window
step first on any new real scenario before a full run, it's much cheaper
than finding out at minute 24 of a full run.

**Comparison result against `hydro.txt` (Cepu gauge, 360h, hourly):**
NSE = -3.83, RMSE = 1432 m^3/s overall -- but that overall number is
misleading because the error is *not* uniform over time:

- **Hours 0-240 (the rising limb): good agreement.** Same shape, ref
  peaks ~2114 m^3/s @ 180h, sim peaks ~2062 m^3/s @ a similar time, errors
  mostly in the 3-8% range. This is the actual validation signal and it's
  a real success -- rainfall-to-runoff generation and routing timing over
  the first ten days are basically right.
- **Hours ~250-360: sim diverges upward and does not turn over.** By
  346h sim is at 5317 m^3/s vs ref's 1966 m^3/s (2.7x). This is *not* a
  NaN-style blowup (everything stays finite, `torch.isfinite` all true)
  -- it looks like a slow-building backwater/storage error that
  compounds over the back half of the event.
- **Only Cepu shows this.** The other 5 gauges (Napel, Kajangan, Sekayu,
  Jurug, Serenan) all peak well before the end of the run and decline
  normally afterward (peaks at steps 1287, 1247, 370, 736, 730 out of
  2160) -- i.e. their hydrographs look physically sane throughout. Cepu's
  peak is at the very last step (2159/2160), i.e. still climbing when the
  run ends.
- **Lead (not yet confirmed as root cause):** Cepu is not the terminal
  outlet -- there are 147 more river cells downstream of it to the real
  sink (river index 139 at grid (54, 291), which also happens to be the
  single widest/deepest channel cell in the whole domain: 150m wide,
  6.6m deep). The bed-elevation profile of the last few cells of that
  reach is unusual: flat around zb=3.0-3.12m for a stretch, then a sudden
  ~1m drop to zb=2.0m in the last two cells -- a real slope estimate over
  the last 5 cells comes out to ~2.2e-4, about **4.5x flatter** than the
  `RiverParams.outlet_slope` default (1e-3) that `river_rhs`'s free-outflow
  boundary formula assumes at the sink cell. Whether/how much this
  mismatch (or something else in that terminal reach) is *causing* the
  backwater-like growth at Cepu was not confirmed before this session's
  time budget ran out -- **this is the next thing to check**, e.g. by
  rerunning with `outlet_slope` set to the real local estimate (~2.2e-4)
  and seeing if the Cepu divergence goes away or changes character, and/or
  recording `qr`/`hr` along the last ~10 cells of the chain (via
  `outlet_riv_index=<list including those indices>`, which
  `RRIModel.simulate` already supports for free via fancy indexing) to
  see whether the growth originates right at the sink and propagates
  upstream, or somewhere else in the reach.

**On the native Fortran reference run (for completeness, not on the
critical path anymore):** compiled the real CUI Fortran source natively
on Linux with gfortran (see below) to get an independent, from-scratch
`out/hs_*.out`/`qr_*.out` reference -- but its `rain.dat` reading step
(a naive double-pass list-directed `READ(*,*)` over ~1.1M tokens, see
`RRI.f90` lines 500-527) turned out to be pathologically slow under
gfortran's list-directed I/O implementation (a known class of gfortran
performance issue, not a bug in our port) -- it was still stuck there
after 45+ minutes of pure CPU-bound parsing (confirmed via `strace`
showing zero syscalls, i.e. it's not stuck on disk I/O, just slow
in-memory token parsing) when this session's validation moved on to using
`solo30s/hydro.txt` instead (ICHARM's own pre-shipped official-run output,
which needs no Fortran re-run at all and was sufficient for the Cepu
comparison above). If picking this up again and a Fortran re-run is still
wanted (e.g. for the spatial `out/hs_*.out` comparison from section 5
step 5, which `hydro.txt` can't provide since it's only one point), either
just let it run in the background for an hour+ and check back, or convert
`rain.dat` to a format read via a tight `do`-loop of individual scalar
`READ`s instead of the implied-do list-directed read (likely much faster
under gfortran), or install ifort/oneAPI if available on the GPU machine
(the original binary's compiler, unlikely to have this issue).

**Building the native Fortran binary (if needed again):** the CUI Fortran
source at `RRI_1_4_2_7_GUI_Beta/RRI-CUI/source/1.4.2.7/*.f90` compiles
clean on Ubuntu's gfortran 11.4 with `-ffree-line-length-none` (ifort
tolerates longer lines than gfortran's free-form default) -- compile all
`*.f90` **except `RRI_Break.f90`** (unused dead code, calls ifort-only
non-standard intrinsics `random_`/`seed_` that don't exist in
libgfortran, and isn't even in the official `make_1_4_2_7.bat` build
list). Two `RRI_Input.txt` values need patching for gfortran's stricter
list-directed integer parsing (ifort accepts a decimal point when reading
an `integer` variable as a non-standard extension; gfortran doesn't):
`riv_thresh` (`100.0`->`100`) and `height_limit_param` (`20.0`->`20`).
Run with `GFORTRAN_UNBUFFERED_ALL=1` in the environment to get real-time
progress out of its `write(*,*) t, "/", maxt` progress line when stdout
is redirected to a log file (gfortran's own I/O buffering is separate
from libc's, so `stdbuf -oL` does *not* work on it, unlike most
C-compiled tools). Also set `OMP_NUM_THREADS=<nproc>` since it's
compiled with `-fopenmp` (`RRI_Slope.f90` and `RRI_Riv.f90` each have a
couple of `!$omp parallel do` regions) -- though note the rain-reading
bottleneck above is entirely serial and won't be sped up by this.

## 6. Known gotchas to remember

- `install/libiomp5md.dll` (already in this repo) needs to be copied next to any
  RRI-CUI `.exe` you run — the binaries are dynamically linked against Intel's
  OpenMP runtime and won't find it otherwise (`error while loading shared
  libraries: libiomp5md.dll`).
- The Read tool misidentifies `obs/*.data` files as binary (they're plain text);
  use `Bash`/`cat` to read them instead.
- `riv_mod.txt` in the Solo project is *not* read by the RRI solver itself (it's
  an intermediate artifact from the GUI's river-extraction tool, `makeRiver3`) —
  don't spend time reverse-engineering its exact meaning; derive the river mask
  from `acc_mod.txt > riv_thresh` instead, per §3.2's table.
- `torch.sqrt`/any `sqrt` of a value that can be exactly 0 needs
  `tensor_ops.safe_sqrt` (floors the argument at a small epsilon), not a bare
  `torch.sqrt(torch.clamp(x, min=0))` — the latter is forward-safe but not
  backward-safe (0 * inf = NaN through unused `torch.where` branches). Already
  fixed everywhere in `rri_torch`; keep this in mind if adding new flux formulas
  for the real-data work (e.g. anything touching `sec_map`/cross-sections).
- Fixed-step RK4 can silently diverge to NaN if `n_substeps` is too small for a
  given scenario's stiffness — this is expected (see README), not a bug to chase.
  Increase substeps before suspecting the physics.
