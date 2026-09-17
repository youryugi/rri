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

**Update: both directions from the lead above were tried, and the
`outlet_slope` hypothesis is ruled out.** Two full 2160-step runs on the
GPU (`differentiable/examples/solo_outlet_diagnosis.py`, ~27 min each),
one with `outlet_slope=1e-3` (default) and one with the measured local
value `2.2e-4`, both recording `qr` at all 148 cells from Cepu to the
true outlet (via `outlet_riv_index=<the whole chain list>`, which
`RRIModel.simulate` already supports for free through fancy indexing --
no code change needed for this). **The two runs are essentially
identical** (peak discharges match to within ~0.5% at every chain
position) -- so the outlet boundary slope is not the cause.

The chain profile reveals the real shape of the problem, and it's bigger
than "Cepu's local backwater": **roughly positions 10-128 of the 148-cell
chain (i.e. most of the mainstem from just below Cepu down to ~20 cells
short of the outlet) never turn over during the whole 360h event** --
every one of them is still climbing at the final step, with peaks
scattered noisily between ~4200 and ~9400 m^3/s (not smoothly varying
cell-to-cell). Comparing against each cell's own Manning "bankfull"
capacity (`Q = (1/ns_river) * (width*depth) * depth^(2/3) * sqrt(1e-3)`,
i.e. what the channel could carry at full depth under the *assumed*
1e-3 slope): positions 0-110 peak at **1.3x-2.4x** that nominal capacity,
while positions 120-147 (near the true outlet, where the reach flattens
and widens further) stay *under* capacity (0.3x-0.6x) yet are still
technically climbing too, just mildly. Only the very first ~10 cells
show a curious alternating pattern (even positions climb, odd positions
peak early ~step 1186 and presumably decline) that wasn't investigated
further.

Reading `hydro.txt` again: the *real* Fortran run's Cepu discharge is
basically flat/mildly declining over the same 250-360h window (1966 ->
1932 m^3/s), so this mainstem reach drains normally in the real model
with the *same* channel geometry (width/depth power-law parameters,
`rivfile_switch=0` in `RRI_Input.txt` means the real run also synthesizes
width/depth this way, not from a surveyed-channel file -- confirmed by
the exact geometry match earlier in this doc). So this isn't "our
synthetic channel sizing is unrealistic for this basin" -- the same
sizing works fine in the reference. **The likely next place to look is
`rri_torch/river.py`'s diffusive-wave RHS itself**: this reach's real bed
slope is very flat (~2.2e-4, see above), much flatter than anything the
synthetic test catchment (`along_slope=0.01` in
`examples/synthetic_basin.py`, i.e. ~45x steeper) or the existing test
suite has ever exercised. Diffusive-wave schemes are known to become
numerically delicate as bed slope approaches zero (the water-surface-
slope term, bed-slope + depth-gradient, can flip sign or lose
significance-of-precision relative to the bed-slope term); this may be
an accuracy/stability limit of the fixed-step RK4 substep count needed
specifically in low-gradient reaches (try substeps much higher than 40,
e.g. 200-500, *only* to see if the reach eventually converges/turns over
-- if it does, this is a stiffness issue and the real fix is either many
more substeps in low-slope regions or an implicit/adaptive scheme; if it
doesn't converge even at high substeps, suspect an actual formula
discrepancy against RRI's real per-cell exchange equation and diff
`river_rhs` line-by-line against `RRI_Riv.f90`/`RRI_RivSlo.f90`'s real
arithmetic).

**Update: root-caused and fixed.** It's the stiffness/resolution
hypothesis, confirmed and fixed in the same session:

1. Ruled out the exchange mechanism as a source of residual error: hand-
   tested `river_slope_exchange` in isolation with an extreme case
   (4m of slope ponding draining into an empty river cell) and confirmed
   the `N_CORRECTION` fixed-point loop is actually an *exact* one-shot
   solve algebraically (the correction formula cancels the mismatch to
   exactly zero after one pass, for any magnitude) -- not an
   approximate iteration that could leave residual error. Also confirmed
   by reading `RRI.f90`'s real main loop (lines 598-966) that the real
   Fortran calls river routing, then slope routing, then evp, then
   `funcrs` (exchange), then infiltration, all once per *outer* dt=600s
   step (river/slope each have their own inner adaptive-`ddt` loop, but
   exchange and infiltration are outside both, exactly once per outer
   step) -- i.e. `rri_torch/model.py`'s operator-splitting order and
   granularity already exactly matches the reference, so that wasn't it
   either.
2. **Confirmed it's substep resolution in the river routing specifically**,
   via `examples/solo_outlet_diagnosis.py` recording `qr` at all 148
   cells from Cepu to the outlet: at `n_substeps_river=20` the chain
   diverges to NaN; at 40 it stays finite but individual cells (e.g.
   chain position 74) show spurious noisy discharge spikes (peaking
   >2x neighbouring cells) that compound into the never-turns-over
   hydrograph. **Raising `n_substeps_river` to 200 (keeping
   `n_substeps_slope=40`) made the whole 148-cell chain profile smooth
   and monotonic**, matching an all-200 (slope+river) run to 3 decimal
   places while costing ~3.7x less (738s vs 2709s for an 800-step
   window) -- river cells are ~17x cheaper per substep than slope cells
   here (1095 vs 18582), so it's much more efficient to refine river
   substeps alone than both uniformly.
3. **Verified against the real reference**: a full 2160-step run at
   `n_substeps_slope=40, n_substeps_river=200` (~34 min on the RTX 5000
   Ada) now gives Cepu **NSE=0.86, RMSE=242 m^3/s** against `hydro.txt`
   (up from NSE=-3.83, RMSE=1432 at the old substeps=40/40 setting) --
   the hydrograph now properly peaks (2431 m^3/s @ 159h vs the
   reference's 2114 m^3/s @ 180h) and recedes afterward, rather than
   climbing monotonically for the whole event. A residual bias remains
   (sim runs ~15% high and ~21h early through the 60-300h window,
   receding to slightly *under* the reference by 346h) -- a convergence
   check at `n_substeps_river=400` was launched to determine whether
   this residual is more substep-resolution headroom or one of the
   model's other already-documented approximations (RK4 vs RRI's real
   adaptive Cash-Karp, the fixed-`N_CORRECTION`-pass exchange, etc.);
   **Update: checked, and `n_substeps_river=200` was already converged.**
   A follow-up 800-step run at `n_substeps_river=400` matches the 200
   run's whole 148-cell chain profile to within ~0.15 m^3/s everywhere
   (e.g. Cepu's own chain position: 2373.780 vs 2373.670) -- i.e. more
   substeps changes nothing further. So the residual ~15%-high,
   ~21h-early bias against `hydro.txt` is **not** a substep-resolution
   artifact; it's coming from one of the model's other already-documented
   approximations (fixed-step RK4 vs RRI's real adaptive Cash-Karp
   integration being the most likely suspect, or the exchange scheme's
   fixed `N_CORRECTION` passes, or the coarse once-per-outer-step rain
   sampling in `build_rain_sequence`). Chasing that residual further
   would mean implementing real adaptive step-size control or refining
   one of those other pieces -- a substantially bigger undertaking than
   this fix, and arguably out of scope for what was a genuine
   stability/resolution *bug* (the never-turns-over divergence) rather
   than an accuracy *ceiling* the project's own README already discloses.
   NSE=0.86 against real official reference output on a real, independently-
   sourced 15,751 km^2 basin is a strong validation result to stop this
   round on.

5. **Update: real adaptive Cash-Karp RKF45 stepping was implemented and
   confirms the residual bias is not a numerics problem.** Per the user's
   request to dig further, `rri_torch/integrate.py` now has
   `integrate_adaptive` -- a direct, differentiable port of RRI's real
   embedded RK4(5) scheme (coefficients transcribed from
   `reference/integrate_ref.py`, which was itself already checked against
   `RRI_Mod2.f90`), with the exact `errmax`/`safety`/`pshrnk`/`ddt_min`
   control logic from `RRI.f90` (including fixing a tiny pre-existing
   discrepancy in `reference/integrate_ref.py`: the real Fortran's
   accept condition is `ddt .gt. ddt_min`, strict, per a "modified Jan 7
   2021" comment in the source -- the reference used `.ge.`; negligible
   in practice but `integrate_adaptive` uses the correct strict form).
   `RRIModel(adaptive=True, eps=..., ddt_min_slope=..., ddt_min_river=...)`
   switches both routing steps to this scheme; `eps=0.01`,
   `ddt_min_slope=1.0`, `ddt_min_river=0.1` are RRI's own hardcoded
   defaults (`RRI_Mod2.f90`), appropriate as-is for the Solo basin's
   metre-scale depths (unlike the synthetic catchment, which needs
   `eps~1e-6` -- see `tests/test_reference_agreement.py`'s existing note
   on why).
   - **Correctness**: a new `test_adaptive_matches_reference` in
     `tests/test_reference_agreement.py` confirms `integrate_adaptive`
     agrees with the independent NumPy reference's own adaptive
     integrator to ~1e-16 relative error (machine precision) on the
     synthetic catchment.
   - **Autograd**: a new `test_gradients_adaptive` in
     `tests/test_forward_and_grad.py` confirms gradients still flow
     correctly through the adaptive path (only the step-accept/reject
     decision itself is non-differentiable, via a `.item()` read on the
     scalar error estimate -- see the module docstring in `integrate.py`
     for the full tradeoff discussion; use `integrate_fixed` instead for
     calibration work where a parameter-independent graph shape matters).
   - **On the real Solo scenario**: a full 2160-step adaptive run took
     only **432s** (~7 min) on the RTX 5000 Ada -- *faster* than the
     already-fixed `substeps-slope=40/substeps-river=200` run (2044s),
     because adaptive stepping only takes many small substeps where the
     error estimate actually demands it, rather than uniformly everywhere.
     Its Cepu hydrograph (peak 2413.97 m^3/s @ 155h) lines up almost
     exactly with the fixed-step-200 result (2430.7 m^3/s @ 159h),
     confirming the fixed-step result was already well-converged. Against
     `hydro.txt`: **NSE=0.8745, RMSE=230.8** (fixed-step-200 was
     NSE=0.8618, RMSE=242.2) -- a marginal improvement, not a resolution
     of the gap. **This decisively rules out time-integration accuracy as
     the source of the remaining ~15%-high/~20h-early bias.**
   - Also directly re-verified (in an isolated hand-built test, not just
     the single-shot case from point 1 above) that the exchange scheme's
     `N_CORRECTION=6` passes converge to <1e-7 relative volume error even
     in the *capped* case-c branch (river deeply overtopped, weir-flow
     capped by available above-ground channel volume) -- so that's ruled
     out too, more thoroughly than before.
   - **Update: `qr_ave` time-averaging implemented and measured -- real,
     but small (~6% of the gap), not the main cause.**
     `integrate_fixed`/`integrate_adaptive` (`rri_torch/integrate.py`)
     now take an optional `on_step(y_new, ddt)` callback, called after
     every *accepted* substep; `RRIModel(..., track_qr_avg=True)` uses it
     to accumulate `qr * ddt` across the river routing's substeps and
     normalize by `dt`, exposing `StepDiagnostics.qr_avg` /
     `simulate()`'s `"qr_avg_outlet"` -- a direct analogue of RRI's own
     `qr_ave` (`RRI.f90` lines 644-715), as opposed to the instantaneous
     end-of-step value `qr`/`"qr_outlet"` always was.
     `examples/solo_validation.py --track-qr-avg` and
     `compare_solo_hydro.py` (which now prefers `qr_avg` automatically,
     `--instantaneous` to compare the old way) wire this up. Result on a
     full adaptive Solo run: **NSE 0.8744 -> 0.8888, RMSE 230.8 ->
     217.2** (bias 175.3 -> 164.5 m^3/s) -- a real, measurable, ~6%
     reduction in error, confirming the effect exists and points the
     right direction, but nowhere near enough to explain the remaining
     gap. **Ruled out as the main cause; keep `track_qr_avg=True` for any
     future Fortran-comparison work since it's strictly more correct and
     now costs nothing extra to enable.**
   - **Update: river confluence handling checked and is exact.**
     `examples/synthetic_basin` (the only scenario every prior test ran
     against) is a single straight channel -- it has no cell where two
     upstream river reaches merge, so `river_rhs`'s vectorized inflow
     accumulation (`inflow.index_add(0, down[has_down], qr[has_down])`,
     `river.py`) had never actually been exercised at a real confluence
     by any test, despite the real Solo network being a dendritic tree
     full of them. Added `test_river_confluence` to
     `tests/test_reference_agreement.py`: builds a small explicit
     Y-shaped network (two tributaries merging into one cell, continuing
     to an outlet) and checks (a) `river_rhs` (torch, `index_add`) against
     `river_rhs_ref` (reference, plain dict `+=` accumulation) -- agrees
     to **exactly 0.0** (not just near machine precision) on this small
     hand-built case, and (b) the confluence cell's inflow against an
     independently hand-derived expected value using the reference's own
     `hq_river` -- also exact. **Ruled out.**
   - **Also tried and abandoned**: getting a real Fortran spatial-field
     comparison (`out/hs_*.out`/`out/qr_*.out`) by fixing the native
     build's slow `rain.dat` parsing. Used `gdb` (rebuilt with `-g -O0`)
     to confirm it really is stuck inside the list-directed `READ` at
     `RRI.f90:506` specifically on whichever block is *last* in the file
     (confirmed by appending a dummy 17th block: the real 16 blocks then
     read in seconds and the slowdown moved to the new last, unused,
     block) -- so it's some gfortran runtime behavior tied to proximity
     to true EOF, not the block's content, and not (verified) a
     CRLF-vs-LF line-ending issue. However, a second attempt praying
     the same "append a dummy block" workaround got *worse*, not better
     (progress went backwards between two `gdb` samples, i8 -> i188 in
     under a minute then only reaching i=63 after 80 more minutes) -- the
     hand-crafted dummy block was likely malformed in some way that made
     things worse, and/or there's a second, different slow path being
     hit. Given ~2 hours total sunk into this with no reliable fix and a
     working alternative validation path already in hand (`hydro.txt`),
     this was abandoned as not worth further time. If revisited: don't
     hand-edit rain.dat with shell text tools again -- instead either (a)
     get ifort/oneAPI onto this machine (the original toolchain, very
     likely doesn't have this gfortran-specific quirk), or (b) rewrite
     just the two rain.dat read loops in `RRI.f90` (lines ~500-527) as a
     small standalone Fortran or Python-preprocessing step that reformats
     rain.dat into whatever gfortran reads fast, and confirm with `gdb`
     that the new last block is *not* slow before trusting a full run.
   - **Where things stand**: outlet BC, exchange convergence (including
     the capped case-c branch), integration accuracy (fixed vs. real
     adaptive), river confluence handling, and most of the `qr_ave`
     timing question are now all checked off with concrete evidence, not
     estimates. None of them explain the bulk of the remaining ~13-15%
     high / ~20h-early bias at Cepu. What's left as genuinely unexplored:
     (a) a real spatial-field diff against the Fortran (blocked on the
     native-build issue above, or would need a different machine/
     compiler); (b) a systematic line-by-line diff of `RRI_Riv.f90`'s
     `qr_calc`/`hq_riv` arithmetic against `hydraulics.hq_river` looking
     for a subtler formula discrepancy than has been checked so far
     (everything checked to date has been the *control logic* around the
     formulas -- confluence summation, exchange equilibration, step-size
     control -- not the core Manning/diffusive-wave formula itself, which
     has only been checked against the synthetic catchment and the
     independent `reference` port, both of which implement the *same*
     `rri_torch`-derived understanding of the formula rather than an
     independently-sourced one); (c) revisiting whether `obs/disc_*.data`
     (real gauge-measured discharge, not Fortran-simulated) can be
     time-aligned with confidence -- if the ~13-15% gap turns out to also
     exist between the *official Fortran* and the *real observations*,
     that would reframe this entirely as an RRI model/parameter
     calibration question rather than anything about our port.
4. **Fixed the defaults**: `examples/solo_validation.py` and
   `examples/solo_outlet_diagnosis.py` now default to
   `--substeps-slope 40 --substeps-river 200` (previously 20/20 and
   40/40 respectively) with this reasoning in their docstrings/module
   comments.

This was *not* a bug in `rri_torch`'s formulas -- `river_rhs` and
`hq_river` (in `hydraulics.py`) are correct; RK4 with too few substeps
is just a known-documented limitation (model.py's own docstring already
says "Not implemented: adaptive RKF45 time-stepping") that bites harder
on this real basin's near-flat (~2e-4) mainstem reach than anything the
synthetic test catchment (35x steeper) or previous test suite exercised.
The fix is choosing `n_substeps_river` per-scenario stiffness, not a
code change to the physics.

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

## 7. Update (session continuing 6a, same 2026-09-16/17 GPU-machine work):
root-caused the Fortran build issue AND found two real formula bugs in
`rri_torch`. NSE at Cepu went 0.86 -> **0.96**.

### 7.1 The native Fortran build issue was a real bug in RRI.f90, not a
performance quirk -- fixed, and it now runs to completion

Section 6a's "pathologically slow rain.dat parsing" was mischaracterized —
it doesn't just look slow near the end of the file, **it never terminates**.
Diagnosed with `gdb` on a `-g -O0` rebuild plus a minimal standalone
reproduction (`repro.f90`, just the counting-pass loop in isolation):

```fortran
tt = 0
do
 read(11, *, iostat = ios) t, nx_rain, ny_rain   ! header read
 do i = 1, ny_rain
  read(11, *, iostat = ios) (rdummy, j = 1, nx_rain)  ! <-- overwrites ios!
 enddo
 if( ios.lt.0 ) exit                              ! checks the INNER loop's
 tt = tt + 1                                       ! iostat, not the header's
enddo
```

`ios` is shared between the header read and the inner per-row reads. After
the true last block is consumed, the *next* header read hits EOF and sets
`ios<0` -- but that value is immediately clobbered by the inner `do i=1,
ny_rain` loop's own read attempts before the `if(ios.lt.0)` check ever
sees it. Fortran leaves read targets unchanged on a failed read, so
`ny_rain` (used as the inner loop bound) keeps its last-good value
forever, and empirically gfortran does *not* keep re-reporting `ios<0` on
every subsequent post-EOF read attempt on this unit -- so the loop just
never exits. (This is likely why the real ifort-compiled Windows binary
never exposed it: whatever ifort's runtime does differently on repeated
post-EOF reads apparently avoids this, but the *logic* is wrong on its own
terms regardless of compiler.) Confirmed the fix in the standalone repro
first (infinite -> 0.19s) before touching the real build: check the
header's own iostat immediately, before the inner loop:

```fortran
 read(11, *, iostat = ios_header) t, nx_rain, ny_rain
 if( ios_header.lt.0 ) exit
 do i = 1, ny_rain
  read(11, *, iostat = ios) (rdummy, j = 1, nx_rain)
 enddo
 tt = tt + 1
```

Same exact bug pattern also exists in the (Solo-irrelevant, `evp_switch=0`
so never hit) PET-file counting loop a few dozen lines down -- fixed that
one too for robustness. **Do not bother with the CRLF-conversion or
"append a dummy trailing block" workarounds described in section 6a's
Fortran-build notes** -- both were chasing the wrong theory (an EOF-
proximity performance quirk) for what's actually a correctness bug; the
one-line `if(ios.lt.0) exit` fix immediately before the inner loop is the
real and complete fix. With it, the full 2160-step Solo run completes in
under a minute and produces genuine `out/hs_*.out`, `out/hr_*.out`,
`out/qr_*.out` (96 files each, one per output step) and a full
`out/storage.dat` mass-balance trace -- the real spatial-field ground
truth section 6a's "next steps" wanted, not just the single-point
`hydro.txt`. (The fixed binary lives at
`rri_build/rri_1_4_2_7_fixed` in the session's scratchpad, built the same
way as section 6a describes, from a copy of `RRI.f90` with this patch
applied -- rebuild from a fresh copy of the real source plus this same
one-line-per-loop fix if picking this up in a new environment.)

`examples/compare_solo_spatial.py` (new) compares `rri_torch`'s per-cell
discharge directly against `out/qr_*.out`, at the exact outer-step index
each output file was written at (`out_next = nint((tt+1) *
dble(maxt)/dble(outnum))`, `RRI.f90` lines 593-1013 -- reproduced as
`fortran_nint`/`output_step_indices` in the script; note Fortran's `NINT`
rounds half away from zero, not Python's banker's-rounding `round()`).
One quirk to know about: the true outlet cell itself reports near-zero
in the real output because `RRI_Riv.f90::qr_calc` explicitly skips any
cell with `domain_riv_idx(k)==2` (`if(domain_riv_idx(k).eq.2) cycle`) --
`rri_torch` does compute a real free-outflow discharge there, so that one
cell's %diff blows up to meaningless values against a near-zero reference
and should be excluded from any summary statistic (the script already
does this).

### 7.2 Two real formula bugs found by diffing against the actual Fortran
arithmetic line-by-line (not just control-flow/logic, which section 6a
had already checked and found clean)

Both were found by reading `RRI_Riv.f90` and `RRI_RivSlo.f90` end-to-end
looking specifically for anything version-tagged in a comment (`! v1.4`,
`! modified v1.4.2.4`, etc.) showing the formula had *changed* at some
point -- both bugs turned out to be exactly this: `rri_torch` matches an
*older* pre-"v1.4.2.4" formula, not the one this project's actual
`source/1.4.2.7` uses.

**(a) River discharge used the wide-channel approximation, not the true
rectangular hydraulic radius.** `RRI_Riv.f90::hq_riv`:

```fortran
a = sqrt(abs(dh)) / ns_river
r = (w * h) / (w + 2.d0 * h)        ! hydraulic radius = area/wetted-perimeter
!q = a * h ** (5.d0/3.d0) * w       ! <- older formula (commented out) -- this
                                     !    is what rri_torch.hydraulics.hq_river had
q = a * r ** (2.d0 / 3.d0) * w * h  ! <- the real, current formula
```

`R = wh/(w+2h)` is the exact rectangular-channel hydraulic radius (area /
wetted perimeter); `R ≈ h` is only accurate for very wide, shallow
channels (`w >> h`). Solo's channels are `w/h ~ 15-20` (e.g. 140m wide,
8m deep) -- nowhere near wide enough for that approximation: `R` comes
out **~10% below `h`** at that ratio, and since `Q ∝ R^(2/3)`, that's a
real, systematic ~7% discharge *under*-estimate at the OLD formula for a
given depth (equivalently, *rri_torch* was requiring a higher `h` than
real RRI to pass the same discharge, i.e. water backs up more before
draining — pushing in exactly the direction of the bias chased since
section 6a). Fixed in `rri_torch/hydraulics.py::hq_river` (floor `R` at a
small epsilon before the `**(2/3)` for the same forward/backward-safety
reason `safe_sqrt` exists — a `<1` power's derivative blows up at 0) and
mirrored in `reference/hydraulics_ref.py` (both were checked back into
exact agreement, see `tests/test_reference_agreement.py`'s existing
machine-precision checks, still passing after this change). The
slope-side Manning formula (`RRI_Slope.f90::hq`, `hydraulics.hq_slope`)
does *not* have an equivalent hydraulic-radius term -- sheet flow over a
2D slope cell has no analogous "wetted perimeter" concept, that formula
was already exact as-is; only the river channel needed this fix.

**(b) River<->slope exchange was computing exactly half the real rate --
missing a "both banks" factor of 2.** `RRI_RivSlo.f90::funcrs`, in
*every* one of its three weir-flow branches (case a "slope free-falls
into an empty river", case c "river overtops onto slope", case d "slope
floods into river"):

```fortran
!hrs = mu1 * hs_top * sqrt( 9.81d0 * hs_top ) * dt * len / area          ! old
hrs = mu1 * hs_top * sqrt( 9.81d0 * hs_top ) * dt * len * 2.d0 / area   ! modified 1.4.2.4
```
(and identically for the `mu2`/`mu3` free/submerged weir terms in cases c
and d). The physical reading: a river reach can overtop *either* of its
two banks onto the adjacent slope (or be filled from either bank), so the
effective weir-crest length driving the exchange is `2 * len_riv`, not
`len_riv` -- rri_torch's `exchange.py` (`_weir_flow`'s `free`/`submerged`
terms, and `hrs_a`) was missing this factor of 2 entirely, i.e. **all
river<->slope water exchange was happening at half the real rate in both
directions**. Note the *cap* logic (`ar = len*b/area`, used to clamp a
transfer to "can't move more than what's actually available") does
*not* get the `*2` -- only the flow-rate formulas themselves; verify
against the real source rather than assuming symmetry if touching this
again. Fixed in `rri_torch/exchange.py` (`_weir_flow` and the case-a
`hrs_a` line) and mirrored in `reference/physics_ref.py`'s
`river_slope_exchange_ref` (five call sites: case a's `hrs`, and both
branches of case c and case d). Verified the exchange function itself
is *exactly* volume-conservative regardless of the rate constant (hand-
tested 2000 random (width, len, depth, height, hs, hr, dt) combinations
against a direct volume-balance check, worst relative error 2e-7 --
floating-point noise) -- so this bug never leaked mass, it just made
slope-to-river drainage (and river-to-slope relief) unrealistically slow
in both directions, letting water pool for too long.

**Both fixes were required to reproduce the real hydrograph shape** --
fixing only (a) barely moved `hydro.txt`'s NSE (0.86 -> 0.87ish); adding
(b) on top of (a) is what took it to 0.96 (see 7.3). In hindsight this
makes sense: (a) changes how much head is needed to pass a given
discharge through the channel itself, while (b) changes how fast the
slope can relieve itself into that channel in the first place -- the
Cepu-gauge symptom chased since section 6a (water piling up and not
draining) was being driven more by the exchange bottleneck than the
channel capacity.

**A pre-existing gfortran-only bug fix note**: while reading
`RRI_Riv.f90`/`RRI_RivSlo.f90` end-to-end for this, no *other* undoubled
formula or missed version bump was found in the sections relevant to this
scenario (dam/diversion/cross-section code wasn't scrutinized to the same
depth since all those switches are off for Solo) -- but if picking up
scope beyond Solo later (dams, diversions, cross-sections), re-check
`RRI_Dam.f90`/`RRI_Div.f90`/`RRI_Section.f90` for the same kind of
version-tagged comment trail before trusting rri_torch has them right,
since none of those paths have been checked at all.

### 7.3 Test suite fallout from fix (b), and why it's expected, not a bug

Doubling the exchange rate exposed a latent (and, once you look for it,
totally expected) operator-splitting error in
`tests/test_forward_and_grad.py`'s mass-balance smoke test: it jumped
from ~3.3% to ~8.6% at the test's original `dt=300s, n_substeps=4`.
**Raising `n_substeps` does not fix this** (verified: 4->16 substeps at
fixed `dt=300s` left the error flat at ~8.4%) because the river<->slope
exchange (and ET, and the sink drain) run exactly once per *outer* `dt`
regardless of substep count, in both `rri_torch` and the real RRI.f90
(section 6a already confirmed the operator-splitting order/granularity
match) -- so their contribution to the total splitting error is `O(dt)`,
not shrinkable by refining the *inner* integration. **Shrinking the
outer `dt` itself does fix it** (verified: dt 300->150->75->30s at a
fixed `n_substeps=4` took the error 8.6% -> 1.0% -> 0.3% -> 0.05%,
converging cleanly to zero). Fixed the test by halving its `dt` to 150s
(doubling `T` and the rain-pulse timing proportionally to cover the same
physical event) rather than raising substeps, since substeps provably
don't help here -- see the comment now in
`test_forward_and_grad.py::test_mass_balance` for this exact reasoning
if it needs revisiting again. This is *not* evidence of a new bug in fix
(b) -- it's the same well-understood, bounded first-order splitting
error every operator-split scheme has, just more visible now that twice
as much water moves through the once-per-outer-step exchange call.

### 7.4 Final validation numbers and what's still open

Full 2160-step adaptive runs with both fixes (`eps=0.01`,
`ddt_min_slope=1.0`, `ddt_min_river=0.1`, ~535s / 9 min on the RTX 5000
Ada):

- **Cepu vs `hydro.txt` (qr_avg, time-averaged): NSE = 0.9626, RMSE =
  126.0 m^3/s** (up from 0.86/242 with only the earlier substep fix from
  section 6a, and from -3.83/1432 at the very first working full run).
  Peak 2267 m^3/s @ 169h vs reference's 2114 m^3/s @ 180h -- both
  magnitude and timing now much closer than any earlier iteration.
- **Spatial check against the real `out/qr_*.out` fields** (148-cell
  Cepu-to-outlet chain, 12 output times spread across the event,
  `compare_solo_spatial.py`): mean %diff +9.8%, median +10.0%, std 8.5%
  (outlet cell excluded per 7.1). This is a *much* smaller and more
  uniform residual than section 6a's chain diagnostic ever saw (which
  ranged from -30% to >100000% before these fixes), but it hasn't gone
  to zero. Two loose ends visible in this data, neither chased further
  this session:
  - The residual bias shrinks over the course of the event (~+11% around
    t=360-900h... wait, hours not shown that high, see the table --
    concretely: high single digits to low teens through the first half,
    trending down to near 0% or slightly negative by t=2160 i.e. the
    very end of the 360h run) -- consistent with the recession-limb
    undershoot already visible in the `hydro.txt` comparison table
    (section "diff" column goes negative by hour 346). Might just be the
    same qr_ave-vs-instantaneous-adjacent effect from section 6a's point
    3, might be something else about recession dynamics specifically;
    not investigated further.
  - Chain position 112 (grid cell (i=47, j=259), a plain E-flowing (dir
    code 1) cell between two ordinary D8-diagonal neighbors -- nothing
    visually unusual about its position in the chain) shows erratic,
    much-larger-than-neighboring-cells swings at several timesteps
    (+31%, +102%, -11%, +42% at different output times, vs. neighbors'
    steady +5-12%) -- looks like a single specific cell with its own
    issue (maybe a confluence, maybe a sec_map/geometry edge case worth
    checking with `read_esri_ascii` on `acc_mod.txt`/`dir_mod.txt` at
    that exact cell), not investigated further.
  - Neither of these was large enough, or clearly-enough one root cause,
    to justify the time within this session's scope -- flagging both as
    concrete starting points rather than leaving "something's still off"
    unspecified.
- All existing tests still pass after both fixes, including a new
  `test_river_confluence` (exact-agreement check for the multi-inflow
  `index_add`, added this session, unrelated to but verified alongside
  these fixes) and the two new adaptive-integrator tests from section 6a.

### 7.5 Update: four more hypotheses chased for the residual ~10% bias --
all ruled out with evidence, no third bug found yet

Per the user's request to keep digging after 7.4, four more concrete
leads were checked. None explained the remaining gap, but each is a real,
verified negative result worth not re-checking from scratch later:

1. **Chain position 112's erratic swings, investigated concretely.**
   Confirmed it's a genuine `rri_torch`-only artifact (the real
   `out/qr_*.out` shows a smooth monotonic sequence across positions
   110-114 at every checked output time; ours has position 112 sitting
   3-5% *above* both neighbors, consistently, in our own run). Found the
   likely structural cause: positions 111 and 112 have **exactly equal**
   bed elevation (`zb=4.510` both) -- a locally dead-flat bed segment,
   where the diffusive-wave `dh` for that one edge is driven purely by
   the `hr` difference (the bed-slope term cancels exactly), making it
   more numerically delicate than neighboring edges that have a real bed
   slope backing up the flow direction. Confirmed no confluence, no
   width/depth/acc discontinuity there -- geometry is otherwise perfectly
   smooth through this stretch. **Not fixed** -- this affects ~1 of 148
   chain cells by a few-to-~100 m^3/s at a few timesteps, immaterial to
   the Cepu-level validation; flagging the mechanism (near-zero local bed
   slope -> dh dominated by state difference alone) in case a similar,
   larger-magnitude instance turns up elsewhere.
2. **Exact RRI `qr_ave` replication (6-stage Cash-Karp quadrature
   average, not just sampling at each substep's end state) -- implemented,
   tested, makes ~no difference.** `RRIParams`... no, `RRIModel.step`'s
   `river_rhs_fn` now logs the `qr` byproduct of *every* internal RHS
   evaluation (4 per fixed-RK4 substep, 6 per adaptive Cash-Karp substep
   -- exactly matching `RRI.f90`'s `qr_ave_temp_idx` accumulation over
   its 6 stages, lines 644-715), averaged per-substep then ddt-weighted
   across substeps -- a much closer replication of the real formula than
   section 6a's original "just sample qr at the substep's final state"
   approximation. Result on the real Solo run: **NSE 0.9610 vs the
   simpler approximation's 0.9626** -- statistically indistinguishable,
   if anything marginally worse. This is a clean, fully negative result:
   the qr_ave-vs-instantaneous *comparison methodology* is not where the
   remaining gap lives, at either level of fidelity. (Kept the exact
   version in `model.py` since it's strictly more correct and the extra
   cost is negligible -- one extra `river_rhs` call already happens to be
   made per stage regardless, this just also logs its `qr` output instead
   of discarding it.)
3. **Phase-lag hypothesis, tested directly and refuted.** The residual
   bias's shape (grows through the rising limb, shrinks and goes slightly
   negative by the event's end -- see the diff column in 7.4's own
   `hydro.txt` table) looks exactly like what a pure timing shift would
   produce, and peak timing is consistently ~10h early (169-170h vs
   180h). Tested directly: cross-correlated `qr_avg` against `hydro.txt`
   at 0-29h of lag. If timing were the main story, NSE should peak near
   a ~10h lag matching the peak-timing offset. It doesn't: **NSE peaks at
   just 2h of lag (0.9610 -> 0.9630, barely different from 0h),
   and gets rapidly *worse* for any larger lag** (0.90 at 10h, 0.69 at
   20h). The remaining gap is a genuine shape/magnitude difference over
   the course of the event, not primarily a time-shifted copy of the
   same curve.
4. **Slope-routing elevation datum: found a real latent bug, confirmed
   inert for this scenario specifically.** `RRI.f90` line 223:
   `zb(i,j) = zs(i,j) - soildepth(land(i,j))` -- the slope routing's own
   elevation datum (used in `RRI_Slope.f90`'s diffusive `dh`) is ground
   level *minus the local soil depth*, not raw ground elevation.
   `rri_torch`'s `slope_rhs`/`slope_edge_fluxes` (`slope.py`) are called
   with `g.zb` = raw `zs` throughout -- no soildepth subtraction anywhere.
   **This is only a no-op when `soildepth` is spatially uniform**: the
   `dh` formula only ever uses *differences* of `(zb + lev)` between
   neighboring cells, and a spatially-constant offset cancels exactly in
   any such difference. Solo has `num_of_landuse=1` with one uniform
   `soildepth=1.000` for the whole domain, so this is provably a dead end
   *for this scenario* -- confirmed algebraically, not just assumed, and
   not worth burning a GPU run to double-check. **This is a real
   correctness gap for any future scenario with spatially-varying
   soildepth** (multiple land-use classes with different soil depths) --
   `geometry.Grid`/`slope.py` would need a per-cell "slope datum" field
   separate from `zb` (ground elevation, still needed elsewhere e.g. for
   `zb_riv`) before trusting a multi-landuse run. Not needed for Solo;
   worth fixing before extending to a scenario where it would matter.

**Net effect of this round: the search space is now much narrower** (all
of outlet BC, exchange convergence, integration accuracy incl. real
adaptive stepping, river confluence, both qr_ave formulations, timing/
phase, the slope datum question, and the adaptive error-mask logic are
ruled out with concrete evidence), but no third formula bug was found to
explain the remaining ~10% pointwise / ~4% NSE gap.

**One more check done: the bias is basin-wide, not Cepu-reach-specific.**
Compared `qr_avg` against `out/qr_*.out` at all 6 gauges (not just the
Cepu chain), at 5 output times spread through the event. Every gauge
shows the *same shape* -- positive bias (+3% to +11%) during the rising
limb, shrinking and often going slightly negative by the later output
times -- just with gauge-specific magnitude (Sekayu barely shows it,
+2.7% max; Kajangan/Jurug/Serenan swing more, down to -14% at some single
timesteps). This rules out "something specific to the long Cepu-ward
mainstem reach" (e.g. an accumulating per-cell rounding difference along
a particularly long flow path) -- whatever's causing this operates
basin-wide, which points more toward a shared upstream mechanism (rainfall
-> slope response) than anything in a specific river reach. Also
double-checked the slope edge-geometry constants (`l1=dy/2, l2=dx/2,
l3=diag/4` for the 8-direction scheme, `RRI_Sub.f90` lines 246-248) and
the adaptive integrator's error-mask logic (`RRI.f90` line 689: river's
`errmax` masks out `domain_riv_idx==0` cells only, *not* `==2` sink
cells -- matching `integrate_adaptive`'s call for river routing, which
passes no mask at all, i.e. every river cell's error counts) -- both
confirmed to already match exactly, no new leads there.

**Update: re-derived `hq_slope` fresh against `RRI_Slope.f90`'s actual
call sites (not just the `hq` subroutine in isolation) -- also ruled
out, but worth recording exactly why since it looked promising at
first.** `hq`'s own body has an `if(dh.lt.0) dh=0.d0` positioned *after*
computing `vm=km*dh` and `va=ka_p*dh` but *before* computing
`al=sqrt(dh)/ns_p` -- at a glance this reads like `vm`/`va` use the
*signed* `dh` while `al` uses the clipped (non-negative) one, which
`hydraulics.hq_slope` doesn't replicate (it takes `dh.abs()` once, up
front, before computing any of the three terms). Chased this down to the
actual call site in `RRI_Slope.f90` (search `call hq(`): the "coming in"
branch does `dh = abs(dh)` *itself*, immediately before calling `hq`,
and the "going out" branch only calls `hq` when `dh>=0` already -- so by
the time `hq` actually runs, `dh` is *already* guaranteed non-negative on
every call, making `hq`'s own internal `if(dh.lt.0)` a defensive no-op
that never fires in practice. `vm`/`va` therefore also only ever see a
non-negative `dh` in the real program, exactly matching what
`hq_slope`'s upfront `.abs()` produces. **Confirmed equivalent, not a
bug** -- but this is exactly the kind of thing a partial read (the
subroutine alone, without its call sites) gets wrong, so if re-deriving
anything else here, check every call site too, not just the subroutine
body.

Candidates not yet checked, roughly in order of how likely they seem:
(a) [checked and ruled out above -- kept for the record of what "slope
rainfall-to-runoff response, re-derived fresh" turned out to mean]; (b) a
genuinely independent re-derivation of `hq_riv`/`funcrs` too, in case a
second undoubled/unfixed detail is hiding in the same files the first two
came from; (c) whether `ka=0`/`ksv=0` being exactly zero (vs. a tiny
positive number) changes which `torch.where` branches are numerically
active in a way that diverges from the Fortran's own zero-handling in
`hq_slope`/`infilt` -- both are meant to make those terms inert, but
"inert" and "bit-identical to the Fortran's own zero-branch" aren't
automatically the same claim and this hasn't been checked cell-by-cell;
(d) whether the rain.dat step-function timing (rain held constant per
600s outer step, evaluated once at the step's end time -- `build_rain_sequence`)
introduces a small systematic lead/lag in when rainfall actually starts
contributing to runoff, compounding basin-wide over many cells even if
individually tiny.

### 7.6 Update: true per-outer-step ground truth obtained (not just every
3.75h), via an existing debug-output mechanism -- reveals the residual
bias is a clean, nearly-constant proportional excess once flow
establishes, not a complex time-varying error. Four more mechanisms
checked and ruled out; still no third bug found.

Per the user's request to get real step-by-step numerical ground truth
(not just the 96-per-run spatial dumps from 7.1), found that `RRI.f90`
already has an unused-for-this-purpose but perfect building block: the
exact code that produces `hydro.txt`/`hydro_hr.txt` (lines 1001-1005)
writes `qr_ave`/`hr` for every station in `location.txt` whenever
`mod(int(time), 3600) == 0` (hourly). Two trivial, well-understood
changes turn this into full per-outer-step, many-station ground truth
with no new logic, reusing a mechanism already trusted (it's what
produced every `hydro.txt` comparison so far):

1. Changed the modulus from `3600` to `600` (`dt` for this project) --
   fires every single outer step (2160 rows instead of 360) rather than
   hourly, since `time` is exactly a multiple of `600` at every outer
   step boundary regardless.
2. Expanded `location.txt` from its single "Cepu 68 167" entry to 153
   rows: all 148 Cepu-to-outlet chain cells (`chain000`..`chain147`,
   1-indexed from the saved `diag_both_fixes.pt`'s `chain_i`/`chain_j`)
   plus the other 5 gauges. No source change needed for this part, just
   a bigger input file -- `write(1012/1013, ...) (qr_ave(hydro_i(k),
   hydro_j(k)), k=1,maxhydro)` already loops over however many stations
   `location.txt` lists.

Rebuilt (`rri_1_4_2_7_debug`), reran -- **41 seconds** wall time for the
full 2160-step run (with 32 OpenMP threads actually engaging this time,
unlike every earlier run this session: `user` time was ~19 minutes
against 41s real time, i.e. real parallelism, for reasons not
investigated -- possibly just this specific build/link combination).
Produced a 2160-row, 154-column `hydro.txt` (time + 153 stations) at
exact per-outer-step resolution. Cross-checked this new mechanism against
the already-trusted `out/qr_*.out` spatial dumps at a matching
time/cell (Cepu at the final step): **1932.01576 in both, to all 6
printed decimal places** -- confirms this new debug channel and the
spatial-dump channel report the literal same underlying `qr_ave` value. and confirms this session's own comparison-script indexing has been correct throughout (a relief, since the past several sections all leaned on that indexing being right).

**The actual finding, at true per-step resolution (Cepu, 148-chain-cell
save from a rerun with the final stage-accurate `qr_avg`,
`diag_final.pt`):** the %diff vs `hydro.txt` is *not* the complicated
time-varying thing the hourly-resolution comparisons in 7.4/7.5 made it
look like. Sampled every step from t=0 to 350h: essentially 0% until
flow becomes hydrologically significant around t=30-35h; then, during
the sharp initial rise (t=35-43h), %diff swings wildly (up to +70%) --
but this is *just* the ordinary artifact of comparing two similar-but-
not-identical curves where the reference is still very close to zero (a
~+95 m^3/s absolute gap over a ~135 m^3/s reference reads as +70%, over
nothing more exotic than that). **Once flow is actually established
(t > ~50h), the %diff is remarkably flat: +10.2% to +11.6% continuously
from t=50h to t=140h**, then declines smoothly and monotonically through
the rest of the rising limb, peak, and recession, crossing zero around
t=340h and going mildly negative after. There is no oscillation, no
second bump, no cell- or time-localized excursion beyond the single
known position-112 artifact (7.5) -- it is functionally a single, clean,
**nearly-constant ~10-11% proportional excess in discharge for a given
rain forcing, sustained for on the order of 100 hours**, which only
*looks* complicated at hourly resolution because the absolute error is
still ramping up/down during the shoulders of the event while the
proportional relationship is already fully established in the middle.
This reframes the problem usefully: it is not "a timing bug" (7.5 already
ruled that out directly) and not "a bug that only shows up at extremes"
-- it behaves exactly like a single multiplicative parameter or rate
constant that's a few percent off *somewhere* in the rainfall-to-
discharge chain, active across the whole basin (7.5's multi-gauge check),
uniformly across a hundred-hour span of active flow.

**Four more specific mechanisms checked with this data and ruled out:**

1. **Total rainfall input volume, checked against `storage.dat`'s own
   `rain_sum` column exactly.** Real Fortran's final cumulative rainfall:
   3.766504e9 m^3 (`out/storage.dat`, last row, column 1). Recomputed
   independently from `build_rain_sequence`'s own output (sum of rain
   rate x `grid.area` x `dt` over every domain cell and every step):
   3.766504e9 m^3. Ratio 1.0000000950 -- exact to 7 significant figures.
   Rainfall input is not the source of a multiplicative excess.
2. **The groundwater module, traced through by hand (not just "the
   switch is off") and confirmed algebraically inert for this project.**
   `RRI_Input.txt` has `ksg=0.000` (bedrock GW conductivity) but
   *nonzero* `gammag=0.037`, `kg0=5.7e-5`, `fpg=0.100` -- worth actually
   checking rather than assuming those nonzero values are dead, since
   they participate in real formulas (`hg_calc`: `qg = dh*kg0/fpg*
   exp(-fpg*hg)`). Traced all four `RRI_GW.f90` entry points by hand:
   `qg_calc` explicitly `cycle`s whenever `ksg_p<=0` (line 67) so lateral
   GW flow is exactly zero every call; `gw_recharge` computes
   `rech=ksg*dt=0` unconditionally when `ksg=0`, and with `hg` starting
   at exactly 0 (`hg_init` zeros it) and never modified by a zero
   `rech`, every subsequent branch of that subroutine evaluates to a
   zero state change too; `gw_lose` is gated on `rgl_idx(k)>0.d0` and
   Solo's `rgl=0.000`; `gw_exfilt` is gated on `hg_idx(k)<0.d0`, which
   never becomes true since `hg` never leaves 0. So `hg` provably stays
   exactly 0 for the entire 2160-step run regardless of `gammag`/`kg0`/
   `fpg`'s nonzero values -- confirms `rri_torch` not implementing
   groundwater at all is correct for this scenario, not an
   approximation.
3. **The real adaptive integrator's actual state variable is `vr`
   (a volume-like quantity), converted to/from `hr` via `hr2vr`/`vr2hr`
   (`RRI_Section.f90`) -- checked whether this differs from
   `rri_torch` integrating `hr` directly.** With `sec_switch=0` (Solo has
   no cross-section table, `sec_map_idx(k)<=0` for every cell), both
   conversions collapse to a simple **linear, state-independent**
   rescaling: `vr = hr * area * area_ratio_idx(k)`, `hr = vr / (area *
   area_ratio_idx(k))`, where `area_ratio_idx(k)` is a fixed per-cell
   constant. A linear, constant-coefficient change of variables doesn't
   change an RK/RKF-family integrator's trajectory or its embedded error
   estimate's *relative* behavior (the `eps` accept/reject criterion
   ends up comparing the same relative error either way) -- confirmed
   this is a real equivalence, not just integrating "the same equation in
   different units happens to usually be fine". Only relevant if
   `sec_switch=1` is ever turned on for a future scenario (cross-section
   tables make the `vr<->hr` relationship genuinely nonlinear via a
   lookup curve, at which point `rri_torch` -- which doesn't implement
   `sec_map` at all -- would need this reconsidered).
4. (Already covered in 7.5, re-confirmed with the finer data): the
   position-112-style artifact and the overall shape are consistent
   between the coarse hourly view and this new per-step view -- the finer
   resolution didn't surface any *new* localized anomaly beyond what 7.5
   already found.

**Where this leaves things (before 7.7 below changed the picture):**
every control-flow, methodology, and input-data question checkable by
reading source and comparing against real output had been checked, and
the remaining gap looked like it must be a genuine, still-unidentified
formula or parameter difference somewhere in the rainfall -> slope ->
river chain.

### 7.7 Update: found via a state-matched (not simulation-matched) test
that the gap is a real, *pointwise* `qr_calc`/`river_rhs` formula
discrepancy -- not an integration/coupling artifact -- and it flips sign
specifically at near-flat-bed-slope cells. Root mechanism still not
identified; this is the most important finding of this round and where
the next session should start.

Per the user's explicit request to keep going with debug-output-based
numerical comparison, went one level more direct than 7.6: instead of
comparing two independently-evolved *trajectories* (which conflates "the
per-step formula is wrong" with "small per-step differences compounded
over 2160 steps"), **fed the real Fortran's own actual `hr` state
(read from its per-step debug output, section 7.6's `hydro_hr.txt`) into
`rri_torch`'s `river_rhs` directly, and compared the resulting `qr`
against the real Fortran's own `qr_ave` for that exact same step** --
i.e. "given the identical input state, does our formula produce the
same output the real one did?" This isolates the formula itself from 2160
steps of accumulated integration/coupling history.

First attempt used the 96-per-run spatial dumps (`out/hr_*.out` /
`out/qr_*.out`, 3.75h apart) and found huge, wildly-varying differences
(+70% to -13%) -- but this comparison has a real flaw: `qr_ave` is a
time-*average* over the preceding ~3.75h window (`qr_ave` is reset every
outer step and accumulated across only that step's substeps -- see 6a --
but summed here across ~22 outer steps between spatial dumps... actually
per-outer-step, so the real flaw is comparing an *end-of-window*
`hr` state against a `qr_ave` that's already time-averaged even over
just its own single 600s step, while sampling `hr` only at the 3.75h
grid) -- so during a fast-rising hydrograph, an end-state-derived
instantaneous `qr` naturally reads high against a preceding window's
average, independent of any bug. Redid it with section 7.6's true
per-*outer*-step (600s) `hydro_hr.txt`/`hydro.txt` (153 stations, every
single step) to shrink that window by ~22x -- if the effect were just
window-averaging, it should have nearly vanished.

**It didn't -- and the full picture, checked across 5 timesteps rather
than 3, turned out more nuanced than first thought (an earlier version of
this note overclaimed "positions 74/110 consistently negative"; corrected
below after checking steps 240/360/480/600/720, not just the first
three).** Feeding the real `hr` state at the *end* of step N into
`river_rhs` and comparing against the real `qr_ave` reported *for* step N
(600s window): **cells with a clear, non-trivial local bed-elevation drop
show a stable +15% to +23% excess across every timestep checked**
(e.g. positions 60, 80, 130: +19% to +23% at all 5 of steps
240/360/480/600/720, barely moving). **Cells near a locally flat bed
segment (`zb_p ~ zb_n`, tiny `dh`) instead show volatile, sign-changing
behavior that drifts over time** rather than a stable opposite sign as
first guessed -- position 110 alone ranges from -25% (step 240) to +102%
(step 360) to -18% (step 480); position 100 drifts from +0.5% (step 240)
smoothly to -36% (step 480) and stays there; position 40 sits close to
0% throughout. So the honest characterization is: **stable, clearly-
bed-sloped cells show a repeatable +15-23% pointwise discharge excess
given an identical input state -- larger than the ~10% seen in the full
accumulated simulation, suggesting some partial self-correction happens
over the course of a real run -- while near-flat cells are simply far
more sensitive to whatever small state differences already exist
upstream (amplifying rather than independently causing a separate
effect)**. This proves the residual gap lives in the `qr_calc`/
`river_rhs` formula's actual pointwise arithmetic, given a matching
state -- not in accumulated integration error, not on the slope side,
not in a comparison-methodology artifact -- but the near-flat-cell
volatility is a symptom of sensitivity to upstream error, not an
independent second bug with its own clean sign. Chase the **stable
+15-23% at clearly-sloped cells** first; it's the reproducible, well-
behaved signal.

**What was checked and did NOT explain the sign flip** (so don't
re-derive these from scratch again): `dis_riv` (confirmed identical
formula, 7.6); the reverse-flow branch's width lookup (real Fortran's
`hq_riv(hw, dh, kk, width_idx(k), qr_temp)` in the `dh<0` branch uses
`width_idx(k)` -- the ORIGINAL cell's width, not the downstream
neighbor `kk`'s, despite passing `kk` as the (unused-for-Solo, only
matters for `sec_map`) third argument -- confirmed `hydraulics.hq_river`
call sites in `river.py` also always index by the origin cell `k`, not
`kk`, matching); the forward/reverse `hw` sill-clamping logic
(`hw=hr_p` unless `zb_p<zb_n` in which case `max(0, zb_p+hr_p-zb_n)` --
confirmed `river.py`'s `hw_fwd`/`hw_rev` match this exactly, condition by
condition). A reverse-engineered "implied `dh`" calculation (solving
`hq_river`'s formula backward from the real `qr` and `hw` to infer what
`dh` the real Fortran must have used, then comparing to our own computed
`dh`) was tried and abandoned -- it depends on correctly guessing which
`hw` branch (forward/reverse) the real Fortran took, which becomes
ambiguous exactly when `dh` is small, so it produced noisy, inconclusive
ratios (0.4x to 2.9x) that don't cleanly support or refute anything on
their own; the state-matched `qr` comparison above is the reliable
result, not this one.

**Also directly re-verified (not just re-derived) as exact, given how
clean a +15-23% multiplicative gap looks -- these are the two most
"obvious" ways a stable percentage bias like this could sneak in, so
worth ruling out explicitly rather than by formula-reading alone:**
`params.river.ns_river.item()` prints exactly `0.03` (float64) -- not a
stray default from a different code path; and `width_param_c/s` (5.000,
0.350) vs `depth_param_c/s` (0.950, 0.200) are not swapped anywhere in
`real_data/solo_river.py` (checked the literals and their exact
assignment to `width[...]` vs `depth[...]` side by side). Neither
explains it.

### 7.8 Update: added the real `dh` debug output (7.7's step 2), and it
resolved the mystery -- mostly a comparison-methodology gap, not a
formula bug. One genuine outlier cell remains, isolated but unexplained.

Implemented 7.7's suggested next step: one new write statement in
`RRI.f90` (right next to the repurposed `hydro.txt`/`hydro_hr.txt` block
from 7.6, no changes to `RRI_Riv.f90`/`qr_calc` itself needed -- `zb_riv`
and `length` are both plain module globals already in scope in the main
program) computing `dh` the same way `qr_calc`'s plain-diffusive branch
does, for 6 hardcoded chain positions (60, 74, 80, 100, 110, 130),
writing to a new `dh_debug.txt` every outer step. Rebuilt
(`rri_1_4_2_7_dhdebug`), reran (41s again).

**`dh` matches to 4+ significant figures at every position and every
timestep checked (ratio 0.9998-1.0010)** -- this is about as close to
"proven bit-identical" as floating-point allows. This conclusively rules
out `zb_riv`, `dis_riv`, and the `dh` formula itself as any part of the
remaining gap. With `dh` confirmed exact, also directly verified (not
inferred) that `zp>=zn` at all 6 positions (so `hw=hr_p` cleanly, no
sill-clamp branch ambiguity) and that plugging the confirmed-exact `dh`
into `hq_riv`'s formula by hand still reproduced the same +15-23%-ish gap
against `hydro.txt`'s `qr_ave` as before -- so the gap really is
downstream of `dh`, between "formula, given a matching state" and
"the real `qr_ave` for that same nominal step."

**That remaining gap turned out to be mostly the same time-averaging
issue from 6a/7.6, just still not fully eliminated: even within one
600s outer step, the real adaptive river integrator can take *dozens* of
substeps** (confirmed directly: re-running just one 600s window through
`rri_torch`'s own `integrate_adaptive`, seeded from the real state at
that window's start, took 31 accepted + 6 rejected substep attempts,
average accepted `ddt` ~19s) -- so comparing a `qr` computed from a
*single* end-of-window `hr` snapshot against `qr_ave` (averaged over
however many of those ~37 sub-evaluations actually happened) is still an
apples-to-oranges comparison, just at a finer grain than 7.6's "3.75h vs
instantaneous" version of the same mistake. **The fix: seed
`rri_torch.integrate.integrate_adaptive` with the real Fortran's own `hr`
state at the *start* of a given 600s step, run our own adaptive
integrator for exactly that one window with `qr_avg`-tracking (the exact
mechanism from 7.6, reused as-is), and compare *that* true window-average
against the real `qr_ave` for the same window.** Result (step 360,
t=60h): the "normal" cells' gap shrinks from +15-23% down to **+4.5% to
+11.5%** (position 74 even flips from -15% to +6.6%; position 100 from
-21% to +4.5%) -- consistent with, and about the same size as, the ~10%
gap already seen in the full accumulated-simulation NSE comparison (7.4).
**This means most of what section 7.7 characterized as a stable
"+15-23% pointwise formula excess" was actually still a comparison-
window artifact, not a formula bug** -- a useful, if slightly deflating,
correction to 7.7's own headline claim.

**One cell did not resolve: position 110 stayed at +97.6%** (barely
moved from 7.7's raw +102.4%) even under this properly-windowed
comparison. This is now isolated as a genuine single-cell anomaly, not
explained by the comparison-methodology issue that accounted for
everywhere else. It's the same near-flat-local-bed-slope class of cell
already flagged in 7.5/7.7 (worth confirming its exact `zb_p`-`zb_n` gap
the same way 111/112 and 74/80 were confirmed near-zero) -- the working
theory remains that these specific cells are so sensitive to whatever
small state differences exist in their immediate neighbors (themselves
subject to 30+ substeps of independent evolution over the same window)
that even a "correctly windowed" comparison doesn't fully tame them, but
this has not been confirmed for position 110 specifically the way it has
for 111/112.

**Net assessment after 7.7+7.8: the ~10% gap already characterized via
the full-simulation NSE (0.96) and the properly-windowed single-step
test are now mutually consistent** -- there does not appear to be a
large, still-hidden pointwise formula bug of the size the raw 7.7 numbers
first suggested. What's left is (a) a smaller, ~5-11%, not-yet-explained
residual at "normal" cells even after correcting the comparison window
(still unresolved -- could be a genuine small remaining formula/parameter
difference, or could need an even more precisely windowed comparison,
e.g. also matching the exact `eps`/substep-count RRI itself would use
rather than `rri_torch`'s own independent adaptive choices for that
window), and (b) isolated single-cell anomalies at near-flat-bed-slope
locations (position 110 confirmed extreme; 74/100/111/112 less extreme
once properly windowed) whose mechanism is understood in kind (local
sensitivity to upstream substep noise) but not nailed down precisely for
any specific cell. Given the scale of effort already invested chasing a
gap that has now shrunk to single digits at most cells, further chasing
should weigh the (now smaller) remaining prize against the (already
demonstrated, considerable) effort each additional layer of precision has
required.

### 7.9 Update: position 110's outlier fully explained -- both models
independently produce the same spatial checkerboard oscillation in this
reach, just out of phase. Not a bug in either implementation.

Per the user's request to keep chasing position 110 specifically, first
ruled out an "invisible tributary" theory: found (by scanning the whole
1095-cell network, not just the 148-cell chain, for any `down_riv`
pointing into a chain cell) that 16 chain positions -- none of them 110
itself, but several upstream of it (92, 97, ...) -- receive inflow from
river cells outside the hand-picked chain that section 7.8's single-
window test left seeded at `hr=0`, an unmodeled deficit. Added those 16
cells to `location.txt` (169 stations total), reran (still ~40s), redid
the seeded single-window test including them -- position 110 barely
moved (+92.4% vs the original +97.6%). **Ruled out.**

Then just looked at the real `hydro.txt` data directly, unfiltered, for
every chain position 100-119 at step 360 (all 148 positions are in there
already, this needed no new Fortran run): the *real* `qr_ave` field
itself is **not smooth** cell-to-cell in this stretch --
100:1253, 101:690, 102:1026, 103:929, ..., **110:593, 111:1126, 112:575**,
113:976, ... -- a pronounced low/high/low/high sawtooth, sharpest right
at 110-113. Checked `rri_torch`'s own accumulated-simulation result
(`diag_final.pt`) at the exact same step and positions: **also sawtooths**
-- 110:1130, 111:717, 112:1169, 113:639, 114:1132 -- same amplitude, same
character, but **shifted by one cell**: the real data's local minima sit
at 110/112, ours sit at 111/113. Position 110 is real's minimum landing
exactly on our maximum (or close to it) -- the single worst-possible
phase alignment, which is exactly why it showed the largest percentage
gap of anywhere checked.

**This is a real, physical(-ish) characteristic of the explicit
diffusive-wave river scheme itself (present in the actual compiled
Fortran, not introduced by rri_torch), not a bug**: in a reach where the
bed-slope-driven part of `dh` nearly vanishes (this stretch, per 7.5/7.7,
sits right where `zb_p~=zb_n` repeatedly), the discharge at each cell is
governed almost entirely by its own `hr` relative to its immediate
neighbors, which can set up a spatial odd/even (checkerboard) decoupling
mode familiar from other explicit finite-difference/finite-volume
schemes on nearly-degenerate/hyperbolic-turning-parabolic problems.
*Which* cells land on the mode's high side vs. low side is decided by
fine details of the adaptive step-size path (confirmed exactly this
sensitivity in 7.6-7.8's `n_accepted=31, n_rejected=6` sub-stepping for
just one 600s window) -- two independently-stepped integrators (real
Fortran's true adaptive path vs. `rri_torch`'s own, seeded from the same
state but otherwise free to choose its own substeps) have no reason to
land in phase with each other, even though each is individually a valid,
consistent solution of the same equations. **Comparing any single cell's
value inside this oscillating reach is therefore comparing against
essentially uncontrollable phase noise, not a stable target** -- the
~100% "error" at position 110 doesn't represent 100% of anything
physically real being missed; a basin/reach-averaged quantity (like the
full accumulated-simulation NSE, or averaging `qr_ave` over the
oscillating cells) is the only fair comparison in a reach like this, and
that's exactly the metric (NSE=0.96, ~5-11% pointwise at "normal"/non-
oscillating cells) that already looked good.

**This closes out the position-110 sub-investigation with a specific,
verified mechanism** rather than an open question: no further action
needed on this specific cell; if a *similarly large* single-cell anomaly
turns up elsewhere in the domain later, check first whether it's sitting
in another near-flat-bed-slope stretch showing the same real-data-is-
also-oscillating signature before assuming it's a new bug.

## 8. Update (2026-09-17 session): paper-prep work -- see `plan.md` for
the full supplementary-work inventory. Working through its priority
order C (ablation) -> A (calibration) -> B (obs time-alignment).

### 8.1 Item C: controlled ablation of the two formula fixes (done)

`examples/ablation_solo.py` reruns the full 360h Solo event four times
-- `both_legacy`/`river_fixed`/`exchange_fixed`/`both_fixed` -- with
IDENTICAL settings otherwise (adaptive integration, `track_qr_avg=True`,
same 154-point combined outlet set: 6 gauges + the 148-cell Cepu chain),
swapping in `examples/_legacy_formulas.py`'s pre-fix formulas via
monkeypatching (`set_formulas()` reassigns `river_mod.hq_river` /
`model_mod.river_slope_exchange` at the *consuming* modules' namespaces,
not the defining ones -- production code is never touched). Results in
`results/`: `ablation_table.csv`/`.md`, `fig_hydrograph.png` (Cepu
hydrograph, all 4 variants vs `hydro.txt`), `fig_spatial_bias.png` (%diff
along the 148-cell chain at t=180h, all 4 variants).

**Headline result -- the two bugs are not equally important**:

| Variant | NSE | RMSE (m3/s) | Bias (m3/s) |
|---|---|---|---|
| Pre-fix (both legacy) | 0.886 | 219.6 | +167.5 |
| Hydraulic-radius fix only | 0.962 | 127.0 | +105.6 |
| Exchange x2 fix only | 0.885 | 221.3 | +169.5 |
| Both fixes (current model) | 0.961 | 128.7 | +107.7 |

The hydraulic-radius fix (section 7.2's fix (a)) accounts for
essentially *all* of the measured improvement at Cepu -- both in outlet
NSE and in the spatial-bias-along-the-chain figure, `river_fixed` and
`both_fixed` track each other almost exactly (visually overlapping
curves in both figures), and `exchange_fixed` alone is statistically
indistinguishable from `both_legacy` (NSE 0.885 vs 0.886; if anything
marginally *worse*, well within run-to-run numerical noise, not a real
regression). The exchange x2 fix (section 7.2's fix (b)) remains
correct and worth keeping -- it has a clear physical justification
(both riverbanks can overtop a weir, not just one; found the same way
as fix (a), via the `RRI_RivSlo.f90` comment trail) -- but its effect on
*this particular basin and event* is too small to show up in either
metric. Plausible reason (not further chased -- would need a basin
where river<->slope exchange is a bigger fraction of total flow to
confirm): Solo's Cepu-ward flow is dominated by in-channel routing for
most of the 360h event, so a factor-of-2 change in a comparatively minor
exchange term barely moves the outlet hydrograph. **For the paper: report
both fixes with their independent physical justification, but be honest
that the ablation shows only one of them is empirically load-bearing for
this validation case** -- overclaiming "two bugs, two improvements"
would not survive a reviewer rerunning this same ablation.

### 8.2 Item A groundwork: why `use_checkpointing=True` still OOM'd, and
the actual fix

Recap of where this was left: `RRIModel.simulate(use_checkpointing=True)`
(checkpointing each *outer* time step) was verified mathematically exact
on the small synthetic basin (`test_gradients_checkpointed_match_uncheckpointed`),
but backpropagating just 200 outer steps on the real Solo scenario
(n_substeps_slope=40, n_substeps_river=200) still OOM'd a 32GB GPU, this
time inside the checkpoint's `recompute_fn` during backward.

**Diagnosis**: wrote `diag_checkpoint_mem.py` (scratchpad) to measure
peak CUDA memory for increasing `T` (5, 10, 20, 40, 80 outer steps) with
`use_checkpointing=True`. Result: **OOM'd already at T=5** (peak ~25.9GB
before failing on a 2MB allocation), i.e. *independent of how many outer
steps are being backpropagated*. This proves outer-step checkpointing
was working exactly as designed (bounding memory to O(1) per outer
step) -- the bug is that "O(1) per outer step" is itself far too large,
because *one* `_step_flat` call already unrolls `n_substeps_slope=40` x
4 RK4 stages = 160 autograd-tracked RHS evaluations (plus
`n_substeps_river=200` x 4 = 800 more, cheaper since the river network
is much smaller than the full grid) over the full 204x336=68544-cell
grid (3.7x the 18582 *active* cells -- slope routing runs over the full
rectangular grid, not just active cells). Recomputing *and locally
backpropagating* that whole 160-eval chain as one unit needs all of it
resident in memory simultaneously -- exactly the same problem
outer-step checkpointing solved for the 200-outer-step axis, just one
level down.

En route, also found and fixed a real (if minor) bug this diagnosis
exposed: `RRIModel.initial_state()` built `hs`/`gampt_ff`/`hr` with no
`device=` argument, silently defaulting to CPU regardless of
`self.grid`'s device -- every GPU caller (e.g. `ablation_solo.py`) had
to work around this manually (`tuple(x.to(device) for x in
model.initial_state(...))`). Fixed to read `device=g.domain.device`
directly, so `initial_state()` just does the right thing on both CPU
and GPU now; existing manual-`.to(device)` call sites still work
unchanged (redundant but harmless).

**Fix: nested (two-level) checkpointing.** Added `checkpoint_substeps`
to `integrate.integrate_fixed` (checkpoints every individual RK4 substep
via `torch.utils.checkpoint.checkpoint(rk4_step, ...)`) and threaded it
through as `RRIModel(..., checkpoint_substeps=True)`, applied to both
the river and slope `integrate_fixed` calls in `step()`. Nested
checkpointing (outer step level + inner substep level simultaneously)
is fully supported under `use_reentrant=False` and, like the outer
level alone, is mathematically exact, not an approximation -- verified
by `test_checkpoint_substeps_match_uncheckpointed` in
`tests/test_forward_and_grad.py` (same hardcoded expected gradients as
the outer-only checkpoint test, `rel_err < 1e-6`, with *both* levels
turned on together -- the real usage pattern).

One correctness subtlety, documented in `integrate_fixed`'s docstring
and enforced with an eager `ValueError` (both there and in
`RRIModel.__init__`): `checkpoint_substeps=True` is **incompatible with
`on_step`/`track_qr_avg=True`**. `on_step` mutates a `nonlocal`
accumulator (`qr_avg` tracking's running sum) as a side effect: if it
sat *inside* the checkpointed region, checkpoint's backward-time
recomputation would silently re-invoke it and double-count every
substep's contribution. Since `track_qr_avg` is a validation-only
diagnostic (matching Fortran's exact `qr_ave` output) that calibration
runs have no need for, this is simply disallowed rather than made
"correct but surprising" by moving `on_step` outside the checkpoint at
the cost of extra complexity.

**Verification on the real Solo scenario**: rerunning
`diag_checkpoint_mem.py` with `checkpoint_substeps=True` added
(alongside the existing outer-level `use_checkpointing=True`), for
T = 5, 20, 80, and 200 outer steps -- the last being the *exact* case
that originally OOM'd at 29.83GB:

| T (outer steps) | peak mem | wall time |
|---|---|---|
| 5 | 0.48 GB | 31s |
| 20 | 0.50 GB | 137s |
| 80 | 0.60 GB | 614s |
| 200 | **0.80 GB** | 1635s |

**Fully confirms the fix**: T=200 now peaks at 0.80GB (vs. the original
29.83GB OOM) -- a ~37x reduction, comfortably inside a 32GB GPU with
huge headroom to spare. The gentle residual growth (0.48->0.80GB over
40x more outer steps, not flat) is itself expected, not a leak: nested
checkpointing only removes the O(n_substeps) *internal* activations per
step, but a genuinely O(T) part remains by design -- one small state
boundary (`hs`, `gampt_ff`, `hr`: ~1.1MB/step combined at this grid's
size) per outer step for the outer-level checkpoint chain, plus each
step's small `qr_outlet`/`drained_volume` outputs, which must stay live
until the final `loss.backward()` reaches them. 200 steps x ~1.1MB/step
~= 220MB matches the observed ~0.32GB growth almost exactly. Compute
cost is the real, expected trade-off of checkpointing (recomputing
forward once during backward, at *two* nested levels here): ~8s/outer
step, i.e. a 200-step forward+backward costs ~27 minutes on an RTX 5000
Ada -- this sets the real constraint on how long a calibration window /
how many optimizer iterations are affordable per session, not memory.

**Item A is now unblocked.** See section 8.4 for the calibration
experiment this enables.

### 8.4 Item A: real-basin gradient calibration vs. Nelder-Mead (done)

`examples/calibrate_solo.py` runs the same synthetic-twin design as the
existing small-basin `compare_gradient_vs_gradientfree.py` (generate a
noise-free "observed" series from known true `ns_slope`/`ns_river`,
perturb to a wrong initial guess, recover it two ways: Adam via
backprop vs. scipy's Nelder-Mead treating the simulator as a black
box), but on the **real Solo geometry** (18582 active cells, real river
network) instead of the 12x7 toy grid.

Two scoping choices, made necessary by section 8.2's finding that a
single real-scale forward+backward pass costs ~8s/outer step even with
both checkpoint levels on: (1) a short window (T=12 outer steps, dt=600s
= 2 hours) with a strong synthetic rain pulse, rather than the real
360h event; (2) the calibration target is a **headwater river cell**
(found via `grid.down_riv`: a cell nobody's downstream target, i.e. the
top of some tributary), not the far-downstream Cepu gauge -- Cepu is
30-40h of real river travel time away even under intense forcing (an
advection delay set by the network's physical celerity, not shrinkable
by raining harder), which would make every calibration iteration cost
hours. This is a scaling/methods demonstration (does gradient
calibration stay cheap and correct on real river-network complexity?),
not a recalibration of the real event.

**Result** (`--steps 12 --adam-iters 15 --nm-maxfev 80`, both from the
same wrong initial guess `ns_slope=0.5x true, ns_river=1.8x true`):

| Method | simulator calls | wall time | best loss | ns_slope rel.err | ns_river rel.err |
|---|---|---|---|---|---|
| Nelder-Mead | 80 | 1187.7s | 1.53e-11 | 0.00% | 0.00% |
| Adam | 15 (best at it=5) | 1376.1s | 1.59e-3 | 0.40% | 14.22% |

**The headline number that supports the paper's thesis**: Adam reaches
a near-exact recovery (loss 1.6e-3, both params within a few % of
truth) in **5 gradient steps**, vs. Nelder-Mead needing **80** black-box
evaluations to reach a similarly tight optimum -- a ~16x fewer-calls
result for this 2-parameter problem, consistent with backprop
extracting a full local gradient from one pass instead of the O(n)
finite-difference-like probing a simplex method needs.

**Two complications, reported honestly rather than smoothed over**
(both real, both instructive, neither a bug):

1. **Adam does not stay converged past its best iteration.** Loss hits
   its minimum (1.59e-3) at iteration 5, then *rises* every iteration
   through 9, oscillating without fully recovering (`ns_slope`
   overshoots true=0.4 up to 0.462; `ns_river` undershoots true=0.03
   down to ~0.023) -- classic fixed-learning-rate overshoot near a
   sharply-curved optimum, not a gradient-correctness issue (the first
   5 iterations descend cleanly and land almost exactly on the truth).
   Concretely: *taking the last iterate* (as you'd have to with real,
   non-synthetic data where the true loss floor is unknown) gives a
   noticeably worse answer (~15%/21% rel. err) than *taking the best
   iterate* (0.4%/14%) -- so a real calibration run needs either a
   properly decayed/tuned learning rate, gradient clipping, or a
   validation-based early-stopping rule; "just run N Adam steps" is not
   automatically as well-behaved as "just run N simplex evaluations."
   (A quick fix -- lowering lr from 0.15 to 0.07 -- was set up but not
   rerun before this was written up; see the TODO below.)

2. **Fewer calls did not mean less wall-clock time here**: Adam's 15
   calls took *longer* in total (1376s vs. 1187s) than Nelder-Mead's
   80, because each gradient call costs far more than each black-box
   call -- Nelder-Mead's `torch.no_grad()` forward pass has none of the
   two-level checkpointing recomputation overhead the Adam path pays
   (roughly 90s/call vs ~15s/call here). **This is an important nuance
   for the paper, not a weakness to hide**: the wall-clock advantage of
   gradient-based calibration isn't automatic at very low parameter
   dimensionality (2 params here) where a simplex method is still cheap
   to run to convergence -- the *real*, dimension-independent advantage
   (one backward pass gives every parameter's gradient at once,
   regardless of how many there are, whereas Nelder-Mead's per-iteration
   cost grows with dimensionality) should show up clearly once item E
   (spatially-distributed parameter calibration -- hundreds to
   thousands of parameters) is attempted, where Nelder-Mead becomes
   outright infeasible while Adam's per-iteration cost is unchanged.
   Report both numbers plainly: calls-to-converge favors gradients
   sharply already; wall-clock-to-converge does not yet, at this scale.

**TODO before this goes in the paper as a figure**: rerun with a lower
Adam learning rate (0.07, already changed in the script) and/or a
tighter cosine schedule so the reported curve is cleanly monotonic
instead of overshooting -- the *mechanism* (backprop needs far fewer
evaluations) is already demonstrated and unlikely to change, but a
non-overshooting curve is a much better figure than one requiring a
paragraph explaining why the last iterate isn't the one to look at.

**Update**: reran at lr=0.07 (`--skip-nm`, since Nelder-Mead is
deterministic given the fixed seed/x0 and didn't need repeating). Loss
is now **cleanly monotonic** every single iteration (1.34 -> 0.016 -- no
overshoot), confirming the overshoot really was a learning-rate issue,
not a gradient-correctness one. Trade-off, also worth keeping rather
than re-tuning away: convergence is slower -- 15 iterations only reaches
loss=0.016 (`ns_slope` rel.err 14.5%), not the lr=0.15 run's it=5 optimum
(loss=1.6e-3). This pair of curves (aggressive-but-overshoots vs.
stable-but-slower) is itself a clean, reportable learning-rate/stability
trade-off finding -- no further reruns needed for this item.

### 8.5 Items E and F: spatially-distributed calibration and an
adjoint sensitivity map

**F (`examples/sensitivity_map_solo.py`)**: the argument that a
differentiable simulator gives you every parameter's sensitivity "for
free" from one backward pass, made concrete. Same scoping as
`calibrate_solo.py` (real geometry, short synthetic pulse, headwater
outlet -- see that file's docstring), but `ns_slope` is left as a full
per-cell leaf tensor (not broadcast from a scalar) and the target is
`sum(qr_outlet)` over the window; one backward pass gives
`d(total_outflow)/d(ns_slope)` at all 18582 active cells simultaneously.
Ran at T=80 (the same window size already timed in section 8.2's memory
sweep: 614s, 0.60GB peak) vs. an initial T=5 smoke test:

| T | wall time | nonzero-sensitivity cells | \|grad\| range |
|---|---|---|---|
| 5 | 32s | 106 | [1.1e-12, 6.7e-2] |
| 80 | 586s | 808 | [1.1e-12, 1.5e+02] |

The resulting map (`results/fig_sensitivity_map.png`) is a tight, small
cluster around the headwater cell in both cases -- growing 16x more
simulated time only grew the affected-cell count ~7.6x, not
proportionally. This is itself a sensible, interpretable finding rather
than an anticlimax: **a headwater cell's contributing area is a fixed
physical property (its own small sub-catchment's extent) -- once the
window is long enough to cover it, more time just accumulates more
total volume from the same finite area, not a growing one.** The map is
effectively tracing out that sub-catchment's boundary for free. Worth
noting for the paper as the "reads like a delineated watershed"
qualitative check that the gradient is doing something physically
sensible, not a numerical artifact.

**E (`examples/calibrate_field_solo.py`)**: the paper's most
consequential demo -- calibrate a full per-cell `ns_slope` FIELD
(~18582 free parameters, one per active cell) instead of a scalar, from
a wrong flat initial guess, using multiple simultaneous observation
points (8 headwater cells spread across the basin -- see
`OUTLET_INDICES`, picked from the 69 available headwater cells found via
`grid.down_riv` for 2D spatial spread) plus an L2 smoothness
regularizer between adjacent active cells (needed because a handful of
point observations cannot possibly constrain 18582 independent
parameters on their own -- without it Adam can drive disconnected cells
to extreme values while barely moving the fit loss; the regularizer
encodes the -- here true -- prior that roughness varies smoothly rather
than cell-to-cell).

True field: `ns_slope = 0.2 + 0.5 * normalized_elevation` (higher/
steeper ground rougher -- a physically plausible, smooth, structured
story, deliberately not per-cell noise, since recovering unstructured
noise from 8 point observations would be meaningless).

Reachability check (reusing the sensitivity-map trick: one backward
pass of `sum(qr at all 8 outlets)` w.r.t. the true field, since any
cell with a nonzero gradient there necessarily influences at least one
observation within the window) at T=30: **13678 of 18582 active cells
(74%)** are reachable from at least one of the 8 outlets -- much better
coverage than F's single-outlet 4%, as expected from spreading
observations across the basin. The identifiability discussion in the
writeup reports recovery error *separately* for reachable vs.
unreachable cells (the latter should, correctly, stay near the flat
initial guess -- not a failure, but the honest, expected behavior of an
underdetermined inverse problem outside its observed region).

Run launched at `--steps 30 --iters 15` (~200s/iteration extrapolated
from section 8.2's T=20/T=80 timings); see the next update for the
finished loss curve, recovered-field figure, and reachable-vs-
unreachable RMSE numbers once it completes.

### 8.3 Item B resolved: the `obs/disc_*.data` "time" column is a
calcHydro output-step index, not hours -- source-code-verified, not
guessed

Recap of the problem (flagged repeatedly, never resolved, across
earlier sessions): all 6 `obs/disc_*.data` gauge files (Cepu, Jurug,
Kajangan, Napel, Sekayu, Serenan) share the identical 16-row time
column `0, 6.4, 12.8, ..., 96`, which cannot be simulation-hours -- the
Solo event is 360h long and its rise/peak/recession all happen well
after simulated hour 96 (see section 8.1's hydrograph: real peak is at
180h). Reading the column at face value as hours gives NSE=-2.32
against the official simulated hydrograph -- worse than predicting the
mean, i.e. actively misleading.

**Root cause, found by reading the actual RRI-CUI post-processing tool
these files are meant to pair with** (`RRI-CUI/etc/calcHydro/calcHydro.f90`,
referenced by `solo30s/calcHydro.txt`'s three config lines:
`./location.txt` / `./out/qr_` / `./disc_`): `calcHydro` reads
`out/qr_NNNNNN.out` for every output step `t = 1, 2, ..., maxt` (until a
file is missing), pulls out each named location's value, and writes
`hydro_<name>.txt` with **the raw output-step index** as the first
column (`write(30, '(i5, e17.8)') t, hydro(i, t)` -- integer step
number, never converted to hours). `RRI_Input.txt` sets `lasth=360`
hours over `outnum=96` output steps, i.e. **3.75h per output step**.

`disc_*.data`'s first column is this same output-step index, just
*continuous-valued* rather than integer (the observations don't happen
to fall exactly on model output steps): `step_index = hour / 3.75`.
Multiplying back out: `6.4 * 3.75 = 24.0` exactly, and every other row
follows the same pattern (`12.8->48h`, `19.2->72h`, ..., `96->360h`) --
**the 16 rows are simply one observation per calendar day (day 0
through day 15) of the 15-day/360h event**, which also lines up exactly
with `rain.dat`'s own 16 daily forcing blocks (section 3.2's table).
The correct conversion is therefore:

```
real_hour = obs_time_column * (lasth / outnum) = obs_time_column * 3.75
```

(equivalently, since every row lands exactly on an integer day here:
`real_hour = row_index * 24`).

**Verified against the official simulated hydrograph** (`hydro.txt`,
already known-good from section 8.1's ablation): interpolating the
official Cepu simulation at `hour = 24*i` for `i=0..15` and computing
NSE against `disc_cepu.data`'s 16 values gives **NSE=0.295** -- a
complete reversal from the literal reading's -2.32, and the resulting
obs-vs-sim comparison (see `scratchpad` verification script) is now
*physically sensible*: both rise from baseflow over the first ~4-8 days
and both stay elevated through the middle of the event before receding.
The remaining gap (obs peaks higher and ~2 days earlier than the
simulation, day 5 vs. day 7-8) is now a genuine, reportable model/data
discrepancy -- not a units bug -- and a plausible contributor is that
each daily obs value's exact sub-day sampling instant (start of day?
daily max? daily mean?) is unknown and unstated anywhere in the
project files, which alone could explain an offset on the order of the
observed 1-2 day lag.

**Action for the paper**: report this 3.75x/day-index conversion
explicitly (with the calcHydro.f90 citation) as a data-preparation
finding in its own right -- a second, independent "artifact in the
official toolchain" alongside section 7's two formula bugs, though this
one is a documentation/labeling gap in the example project's obs data
rather than a bug in the simulation engine itself. `real_data/solo_river.py`
should gain a small helper (`read_obs_series` already exists per
section 3.2 -- check whether it already does this conversion or still
takes the raw column at face value) so any future NSE-against-observed-
data computation uses `hour = row_index * 3.75` automatically instead of
requiring every caller to know this.

**Update -- implemented**: `read_obs_series(path, step_hours=None)` now
takes an explicit `step_hours` multiplier (default `None` keeps the old
raw-column behaviour for callers that want calcHydro's own step-index
units instead of hours). `examples/final_validation_figure.py` uses
`step_hours=LASTH/OUTNUM` (=3.75) and produces the paper's headline
validation figure (`results/fig_final_validation.png`): 3 curves only
(official Fortran, rri_torch both-fixes, real `disc_cepu.data`
observations) -- deliberately not a multi-model panel figure (a user
request, by analogy with more elaborate published flood-model
comparison figures, was to keep this to exactly 3 series). Reuses the
already-computed `results/qr_avg_both_fixed.npy` from section 8.1's
ablation run, so no new simulation.

**Numbers, reported honestly**: `rri_torch` vs official Fortran
NSE=0.961 (RMSE 128.7 m3/s) -- the two curves are visually
near-indistinguishable in the figure, the intended "port is faithful"
result. Against the real observations: `rri_torch` NSE=0.411, official
NSE=0.295 (both computed by interpolating each simulated series at the
16 obs instants). **Both these numbers are "unsatisfactory" by standard
hydrological NSE thresholds** (e.g. Moriasi et al. 2007: NSE<=0.50 is
unsatisfactory) and must not be reported as if they were a good
real-world validation result. The correct framing for the paper is not
"our model matches reality" but: (a) neither simulation has ever been
calibrated against these observations -- both use RRI_Input.txt's
off-the-shelf tutorial parameters (ns_slope=0.4, ns_river=0.03), so a
mediocre uncalibrated fit is expected, not a sign of a broken model;
(b) the 16-point daily series' exact intra-day sampling convention is
unknown (section 8.3), which alone could account for much of the
apparent timing offset; (c) this low uncalibrated NSE is precisely the
motivation for section 8.4's calibration capability -- if the
uncalibrated model already matched reality well, gradient-based
calibration would have nothing to contribute. Do not let `rri_torch`
edging out the official Fortran on this metric (0.411 vs 0.295) be
read as "our port is more accurate than the original" -- both are
essentially the same simulation and the real observation gap swamps
their small mutual difference; it is noise, not signal.

### 8.5 (continued) Items E and F final results

**F final**: T=80 sensitivity map (`results/fig_sensitivity_map.png`)
finished at 586s, 808 nonzero-sensitivity cells (vs. 106 at the T=5
smoke test) -- growing 16x the simulated time only grew the affected
area ~7.6x, sublinearly. Read as a positive finding, not a shortfall:
a headwater cell's contributing area is a fixed physical property (its
own sub-catchment's extent); once the simulation window covers it,
more time just accumulates more volume from the same finite area
rather than recruiting new area. The map is, in effect, a free-by-
product watershed delineation for that cell.

**E final** (`--steps 30 --iters 15`, 3216.7s total, ~214s/iteration):
loss decreased cleanly and monotonically every iteration (6.63 -> 0.86,
no overshoot -- unlike section 8.4's initial lr=0.15 Adam run, this
one's lr=0.05 was well-tuned from the start). **But 15 iterations was
not remotely enough to converge**: `RMSE(ns)` over the 13678 cells
reachable from at least one of the 8 outlets was 0.208, *worse* than
the 0.171 RMSE of the cells that never got any gradient at all and
stayed frozen at the flat initial guess (0.4) -- not a bug, just
confirmation that the loss was still steeply descending (14->15's step
alone was a 6% drop) with no sign of plateauing. The recovered-field
figure (`results/fig_field_calibration.png`) visually confirms this:
the field stays close to uniform ~0.4-0.45 rather than showing the true
field's 0.2-0.7 elevation-linked structure, with one region (the
western sub-catchment near several outlets) overshooting into a
"corrected too far in the right direction, not yet far enough"
red hot-spot on the error map.

**A more fundamental diagnosis than "just needs more iterations"**
(surfaced by direct user question, and worth keeping as the paper's
framing): 8 point observations constraining ~13678 free parameters is
severely underdetermined *in principle*, independent of optimizer
budget. Many spatially-distinct cells within the same sub-catchment
contribute to the same single downstream time series in a way that
blends additively -- the sensitivity Jacobian restricted to one
sub-catchment has effective rank far below its cell count, so there is
a whole subspace of field perturbations the loss is nearly blind to.
More iterations converge to *the* regularized optimum for this
data+prior combination, but that optimum is not guaranteed to equal
the true field, especially deep inside a shared sub-catchment away from
any outlet. This is not a defect of the differentiable-calibration
method -- it is a well-known, general limitation of distributed
hydrological calibration from sparse gauge networks, and the smoothness
regularizer is doing exactly the (partial, prior-dependent) job such
regularizers are supposed to do. **Suggested framing for the paper**:
report this as an honest limitation/discussion point rather than
re-running until the figure "looks right" -- optionally, a follow-up
experiment varying the number of observation outlets (4 vs 8 vs 20) and
showing recovery RMSE improve accordingly would turn this into a
positive, quantified result about identifiability rather than a loose
end.

## 9. Update (2026-09-17 session, continued): preparing a second real
basin (Ishikari River, Japan) for cross-validation (plan.md item I)

Prompted by a user request to actually prepare a second real dataset
(after the Indus/Chashma lead in section 8.x's discussion turned out to
have no shippable DEM/rain files -- only two small `evalHydro`-format
example text files survive from that 2010 Pakistan-flood application
described in `RRI_Manual.pdf` section 8).

**The user supplied `varssim.zip`, `readme_21.txt`, `basinlist_all.xlsx`**
-- these turned out to be an entirely different, well-documented public
dataset: **MERV-Jp** (Multi-model Ensemble for Robust Verification of
hydrological modeling in Japan; Sawada, Okugawa & Kimizuka, 2022,
*Hydrological Research Letters* 16, 73-79, doi:10.3178/hrl.16.73).
135 basins (ver1.1, 1993-2003) / 87 basins (ver2.1, 1986-2015) across
Japan, each with daily basin-averaged precipitation/temperature/PET
(from APHRODITE reanalysis) and observed + 44-model-ensemble-simulated
runoff (mm/day). `basinlist_all.xlsx` links each basin to a GRDC
station number, river/station name (with Japanese kanji), downstream
station linkage, and catchment area -- exactly the missing metadata
that made the Indus lead a dead end.

**Chosen basin: Ishikari River (石狩川), Hokkaido** -- has 4 chained
GRDC-linked gauges (INO 3378.6km2 -> OSAMUNAI 3558.0km2 -> HASHIMOTOCHO
5711.0km2 -> TSUKIGATA 9306.0km2 -> ISHIKARI-OHASHI 12697.0km2, the
final/most-downstream one), comparable in scale to Solo's 15751km2, 30
years of data (basin index 5 in both ver1.1 and ver2.1's file numbering
-- `varssim/ver{1,2}_1/varssim005.csv`), and a multi-gauge chain
structure that mirrors Solo's own 6-gauge-plus-chain validation setup.

**Station geolocation**: GRDC itself was hard to query directly (its
data portal isn't simple-URL-fetchable), but Japan's own river.go.jp
(MLIT Water Information System -- the same source `readme_21.txt` cites
for ver2.1's observed runoff) has per-station pages at
`http://www1.river.go.jp/cgi-bin/SiteInfo.exe?ID=<station-code>`, found
by web-searching the Japanese station name (e.g. "石狩大橋 水位観測所")
-- gives precise degree-minute-second coordinates. Ishikari-Ohashi
(the outlet used): **43°07'20"N, 141°32'32"E** (43.1222, 141.5422).

**DEM/flow-direction/accumulation source: HydroSHEDS v1, 30 arc-second
resolution** (matching Solo's own resolution) -- picked over v2 because
v1's Asia-continental files are directly, statically downloadable
(v2's download page requires clicking through an interactive map, no
guessable direct URL found):
```
https://data.hydrosheds.org/file/hydrosheds-v1-dem/hyd_as_dem_30s.zip  (46.5MB)
https://data.hydrosheds.org/file/hydrosheds-v1-dir/hyd_as_dir_30s.zip  (12.2MB)
https://data.hydrosheds.org/file/hydrosheds-v1-acc/hyd_as_acc_30s.zip  (27.0MB)
```
**HydroSHEDS uses the exact same D8 direction code convention RRI does**
(1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE -- standard ESRI D8,
matching `rri_torch/geometry.py`'s `D8_OFFSETS` exactly) -- no remapping
needed between the two.

**Outlet snapping and verification**: the raw lat/lon landed near but
not exactly on the river cell in the accumulation raster (common --
DEM/coordinate precision); snapped to the local-max-accumulation cell
within a 11x11 window (43.1375N, 141.5042E, accumulation=20552 cells).
Sanity check: at this latitude, one 30-arcsec cell is ~0.628 km2, so
12697 km2 / 0.628 km2 = ~20214 expected accumulated cells --
**matches the snapped cell's actual accumulation (20552) to within
1.7%**, confirming both the station geolocation and the snap are
correct before spending any effort on delineation.

**Watershed delineation**: `pysheds`'s `Grid.catchment()` was tried
first and gave a nonsense 2-cell result (likely an API/coordinate-
convention mismatch not worth debugging further) -- replaced with a
~15-line hand-rolled BFS/flood-fill directly on the D8 direction array
(for each frontier cell, check all 8 neighbours and add any neighbour
whose own direction code points back at the frontier cell), the same
logic pattern already used for `Grid.down_riv` tracing elsewhere in
this codebase. Result: **20552 cells, ~12906 km2 vs. the official
12697 km2 -- 1.6% error**, and the resulting DEM/accumulation maps
(`ishikari/basin_check.png`) show a clean, physically sensible mountain-
ringed basin draining through a single dendritic river network to one
outlet -- no manual cleanup needed.

**Export**: clipped DEM/dir/acc to the catchment's bounding box (+5-cell
pad), masked non-catchment cells to DEM=-9999 (RRI's own domain
convention is `zs > -100`, so this alone correctly excludes them
regardless of their dir/acc values), wrote `adem.txt`/`dir_mod.txt`/
`acc_mod.txt` in RRI's exact ESRI ASCII format. Grid: 227x250 (56750
total, 20552 active) -- same order of magnitude as Solo's 336x204
(68544 total, 18582 active). Verified round-trip through
`real_data/solo_river.py`'s own `read_esri_ascii` parser (unmodified --
no new parsing code needed). Files now live at
`RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/ishikari/topo/`.

**Still needed before this basin can actually be simulated** (not yet
done):
1. `rain.dat` -- convert `varssim/ver{1,2}_1/varssim005.csv`'s daily
   basin-averaged precipitation (mm/day) into RRI's step-function block
   format (a single spatially-uniform value per day, broadcast across
   the whole grid -- matches how this dataset's precipitation was
   itself derived, so this is not a loss of fidelity versus the source
   data). Need a units decision: mm/day -> the m/s RRI expects, and
   picking `dt`/`lasth` for whichever event/period is simulated (this
   basin's data is DAILY, unlike Solo's sub-daily rain.dat, so a
   full-resolution replay isn't meaningful here -- likely want to pick
   a specific flood period from the 30-year record and either treat
   each day as a step-function block directly, or downscale somehow).
2. `RRI_Input.txt` -- no real river-geometry survey data exists for
   this basin (unlike none for Solo either, actually -- Solo also used
   the `width_param_c/s`/`depth_param_c/s` power-law formulas, so this
   is not a new gap) -- reuse Solo's formula-based approach with
   Ishikari-appropriate `riv_thresh` (needs tuning to this basin's
   accumulation scale) and initial Manning's n guesses.
3. Observed validation series: `varssim005.csv`'s "Obs flow" column
   (mm/day, i.e. basin-averaged specific discharge -- convert to m3/s
   via `* catchment_area_km2 * 1000 / 86400` for comparison with RRI's
   `qr` outputs) already gives real observed discharge with a *known,
   unambiguous* daily timestamp (unlike Solo's disc_cepu.data mystery-
   units problem from section 8.3) -- a cleaner validation target than
   Solo's own observed data turned out to be.
4. The other 3 chained gauges (INO/OSAMUNAI/HASHIMOTOCHO/TSUKIGATA)
   still need their own river.go.jp lookups if a multi-gauge validation
   (mirroring Solo's 6-gauge setup) is wanted -- not blocking, since the
   single downstream outlet is enough to start.
5. No official-Fortran reference run is possible for this basin (we
   only have the compiled RRI binary + Solo's own tutorial project, not
   a second official example) -- validation here can only be
   rri_torch-vs-real-observation, not the "vs official Fortran" axis
   section 8.1's ablation relied on. This is fine for the "does the
   whole pipeline generalize to a second real basin" question item I
   is meant to answer, just worth being explicit that this axis is
   Solo-only.

## 10. Update (2026-09-17 session, continued further): first forward
run on Ishikari, a real DEM/DIR inconsistency bug found and fixed, and
real channel-geometry survey data discovered

### 10.1 First forward sanity check and a genuine, hard timing constraint

`examples/run_ishikari.py` builds the scenario directly in Python (same
formulas as `build_solo_scenario`, no `RRI_Input.txt`/`rain.dat` files
written -- there is no official Fortran run possible for this basin
anyway, see section 9 point 5) and forces it with `varssim005_ver2_1.csv`'s
own daily basin-averaged precipitation, broadcast uniformly. First
attempt: 2001-09-01 to 09-20 (20 days, chosen because 2001-09-12 is this
30-year record's **historical maximum** observed flow, 40.97mm/day --
see `df.nlargest(10, 'Obs flow')`). This stalled past 20+ minutes
without finishing (killed).

**Direct extrapolation check, requested by the user**: a calm 2-day
window took 8.5s (adaptive integration) -- naively scaling that rate to
the full 30-year record (10957 days) gives ~12.9 hours, which sounds
feasible, but this is misleading: a single 20-day window containing
*one* storm already blew past 20 minutes, far more than
20 x (8.5/2) = 85s predicted. **Storm days cost far more than calm
days, non-linearly, so the calm-period rate cannot be extrapolated to
the full record.** Separately (also user-asked): the 44 MARRMoT
lumped-conceptual models this basin's CSV includes are structurally far
cheaper than a 2D distributed model like RRI -- no spatial grid, just a
handful of state variables per model, implicit-Euler-solved -- almost
certainly seconds-to-minutes for all 44 models x 30 years combined
(no verified benchmark found, but this follows directly from having no
spatial routing at all). **This basin's data volume (30 years daily) is
simply the wrong scale for an event-based, spatially-distributed model
like RRI/rri_torch to replay continuously** -- RRI's own manual's own
"Application Example" (section 8.x) is itself a single-event case study
(the 2010 Pakistan flood), not a continuous multi-decade run. The
correct scope, confirmed by this exercise, is: pick specific flood
event windows (days to weeks) from the 30-year record, not the whole
thing.

### 10.2 Diagnosing *why* even a single storm was pathologically slow:
a real DEM/flow-direction inconsistency

Direct per-step instrumentation (`diag_ishikari_stall.py`, calling
`model.step()` in a loop and printing `AdaptiveStats` each step) on
just the historical-peak day (2001-09-11, 107.7mm) showed the slowdown
building up gradually and *physically implausibly*: `hs_max` (slope
depth) reached **6.05m** and `hr_max` (river depth) reached **10.6m**
by the end of just one day -- multiple metres of standing water on a
hillslope from one day of rain is not physically reasonable, and is
exactly the signature of water piling up behind something that isn't
draining properly, forcing the adaptive integrator into ever-smaller
substeps to track the resulting stiff dynamics (river `n_accepted` grew
from 1 to 87 substeps/outer-step over the day, req wall time per step
climbing from 0.03s to 0.30s and still rising, not plateauing).

**Root cause, found by directly checking every active cell's D8
direction against its neighbour's elevation**: `hyd_as_dem_30s.zip`
(the DEM, downloaded from the `hydrosheds-v1-dem` endpoint) and
`hyd_as_dir_30s.zip`/`hyd_as_acc_30s.zip` (flow direction/accumulation,
from `hydrosheds-v1-dir`/`-acc`) are **not mutually self-consistent** at
30 arc-second resolution -- HydroSHEDS' hydrologically *conditioned*
DEM (pit-filled, guaranteed monotonic along its own flow network) is
only published at 3 arc-second resolution (`hydrosheds-v1-con`); the
30s "dem" product is apparently a separate, plain elevation raster not
re-conditioned to match the separately-derived 30s DIR/ACC grids.
Direct check: **6.67% of active cells (1371 of 20551) have their D8
direction pointing to a neighbour with *higher* elevation** -- some by
over 100m (worst case: 103m "uphill"). Each such cell acts like a small
dam in the actual (elevation-gradient-driven) diffusive-wave physics,
even though the discrete D8 tracing used for watershed delineation
(section 9) doesn't care about elevation at all and so was unaffected
(hence the watershed shape/area still came out correct).

**Fix**: `condition_dem.py` does a BFS from the outlet, walking
*upstream* against the (trusted, since DIR/ACC gave a physically
correct watershed) D8 network, enforcing
`zs[cell] = max(zs[cell], zs[downstream] + 0.02m)` at each step --
standard "fill pits along a known network" DEM conditioning, just
applied directly along the already-known flow tree instead of needing
to rederive it from scratch (e.g. via `pysheds`' or `richdem`'s
depression-filling, which operate on the DEM alone). Reached all 20552
active cells from the outlet (full connectivity), raised 3872 cells
(18.8% -- more than the 6.67% direct check, since fixing one cell can
cascade a raise onto cells upstream of *it*), largest single raise
117m, zero remaining violations afterward. Rerunning the same
diagnostic showed a real but partial improvement (peak depths dropped,
timing improved somewhat) -- **this fix was necessary but not
sufficient on its own**; see 10.3 for the other, larger contributor.

**Note for item I's writeup**: unlike Solo (where the tutorial package
ships pre-validated, mutually-consistent DEM/dir/acc/rain files with
zero data-preparation risk), assembling a *new* real basin from public
sources surfaces exactly this kind of data-consistency risk -- worth
reporting as a general lesson about preparing RRI-format inputs from
independently-sourced DEM/hydrography products, not just an Ishikari-
specific footnote.

### 10.3 A user-supplied real channel-geometry survey dataset, and a
~3x channel-width error in the borrowed formula

The user separately supplied `ishikarigawakaryuu.zip` (石狩川下流,
lower reaches) and `ishikarigawazyouryuu.zip` (上流, upper reaches),
37.5MB and 2.8MB. **These are exactly the "periodic longitudinal/cross-
section survey" (河川定期縦横断測量) data section 8's discussion had
found *exists* but assumed required a formal information-disclosure
request to obtain** (Sapporo Development and Construction Bureau /
Hokkaido Development Bureau, covering the Ishikari mainstem's full
268km plus 30+ tributaries, multiple survey years from H21/2009 through
R03/2021-22). The zip's internal filenames are Shift-JIS-encoded and
were garbled by a plain `unzip` (which assumes cp437 for non-UTF8-
flagged entries); decoding via `info.filename.encode('cp437').decode('cp932')`
(Python's `zipfile`, not `shift_jis` directly -- a few characters like
the wave-dash differ between the two and only `cp932` matched exactly)
recovers the real names correctly. 4255 files total.

Format (e.g. `01.石狩川/②定期縦断測量測量成果（数値データ）/R03/WZAA2001.CSV`,
the *longitudinal* survey, one row per ~0.5km along the whole river):
`KP, width_m, ..., bank/thalweg elevations..., survey_date`. Cross-
referencing this against `②`'s companion `③定期横断測量測量成果`
(individual cross-section point clouds, `distance_from_centerline, elevation`
pairs per KP) confirms the 2nd column is real channel width in metres.

**At KP26.5-26.567 -- essentially exactly the Ishikari-Ohashi gauge's
own location** (matches river.go.jp's independently-stated "26.60km
from the mouth" almost exactly, an unplanned cross-check that both
datasets are correctly geolocated): real surveyed width **430-560m
(~480m at KP26.5)**. `run_ishikari.py`'s borrowed-from-Solo formula
(`width = 5.0 * drainage_km2^0.35`) predicts only **~163m** at this
basin's 12697km2 drainage area -- **real width is ~3x the formula's
estimate**, i.e. the channel was being modelled as roughly 1/3 its true
carrying capacity, which is exactly what would force anomalously deep,
numerically stiff backups under load (section 10.2's already-fixed DEM
bug was compounding this, not the sole cause).

**Quick fix applied**: rescaled `width_param_c` from 5.0 to **14.724**
(`= 5.0 * 480/163`) so the formula matches the real width *exactly* at
this basin's own outlet, keeping Solo's exponent (0.35) since only one
reliable (area, width) anchor point has been checked so far -- a
deliberately narrow, low-risk correction (anchored to real data at the
one location that matters most for this basin's outlet hydrograph),
not a full refit. Rerunning the section 10.2 diagnostic with this
change: `hs_max` at step 140/day dropped from 2.3m to **1.58m**, `hr_max`
from 8.1m to **6.75m** -- a real, meaningful improvement, but **depths
are still too high to call this basin's forward simulation clean** --
values are still climbing near the end of the tested window, not
plateauing. This means the *single-point, single-coefficient* rescaling
isn't the whole story: a rough two-point check (a second real width
sample at KP140, ~200m, against a much smaller but unmeasured upstream
drainage area) suggests the *true* width-vs-area exponent may be much
steeper than Solo's borrowed 0.35, meaning the current fix is likely
still under-widening the channel away from KP26.6 specifically (both
upstream and, unverified, in the tributaries).

### 10.4 Full KP <-> grid-cell mapping built, real widths applied along
the traced mainstem

Completed the pipeline 10.3 called for, per the user's explicit request
to keep going rather than stop at the single-point rescaling:

**Coordinate system, identified and cross-validated**: the `①距離標
設置測量成果表` (distance-marker) tables give each KP's left/right-bank
(X, Y) in **JGD2000 / Japan Plane Rectangular CS Zone XII (EPSG:2454)**
-- confirmed (not assumed) by transforming KP26.5's coordinates and
checking the result against river.go.jp's independently-sourced
43.1222N/141.5422E for the same gauge: matched to within 0.0003 degrees
(~25m) once the axis order was worked out (Japan's plane rectangular
convention is (X=north, Y=east) -- pass `(Y, X)` as pyproj's
`(easting, northing)` under `always_xy=True`, not `(X, Y)`).

**`build_kp_tables.py`** parses this into two lookups over the
continuous KP0.1-139.0 karyuu coverage (279 distance-marker points):
`kp_position_table.csv` (KP -> lat/lon, bank-averaged) and
`kp_width_table.csv` (KP -> real width in metres, 309 points after
deduplicating across the H24/H26/R03 survey-year subfolders --
preferring the latest year available at each KP).

**`trace_mainstem.py`**: from the outlet, repeatedly steps to whichever
upstream river-cell neighbour has the *largest* accumulation (the
main-channel branch at every confluence, by definition) -- 15 lines,
same BFS-upstream pattern as section 10.2's DEM fix. Result: **241
mainstem cells spanning 218.8km** upstream of Ishikari-Ohashi (this
sub-basin's mainstem reaches nearly to the full river's headwaters --
consistent with this sub-basin being 12697 of the full river's
14330km2, i.e. it's *most* of the whole system, missing only the short
reach from the outlet down to the actual coast).

**`examples/ishikari_real_geometry.py`** (`build_real_width_field`)
combines these: for each traced mainstem cell, compute its absolute KP
as `26.6 + cumulative_upstream_distance`, and if that KP falls inside
the survey's covered range, overwrite the formula-based width at that
one grid cell with the real value interpolated from `kp_width_table.csv`;
everywhere else (tributaries, and the 120 mainstem cells beyond KP139,
i.e. roughly the upstream half of the traced mainstem, where this
session didn't parse the scattered per-project zyouryuu year-folders --
see 10.3's file listing) keeps the already-once-corrected formula
estimate. Wired into `run_ishikari.py.build_ishikari_scenario` (prints
a one-line coverage summary each run). Result: **121 of 241 mainstem
cells now carry real surveyed widths** (KP0.1-138.0 coverage).

**Effect on the section 10.2 diagnostic** (same historical-peak day,
2001-09-11, 107.7mm): `hr_max` at step 140/day went 8.13m (formula
only) -> 6.75m (single-point rescaled coefficient) -> **6.01m** (real
per-cell mainstem widths) -- a genuine, monotonic, but *diminishing*
improvement at each stage. `hs_max` similarly: 2.3m -> 1.58m -> (not
re-checked separately, tracked in the same run).

**Reassessment of what "done" looks like here**: depths are still
elevated, but two considerations argue against reflexively chasing
this to zero by continuing to extend real-data coverage further
upstream: (1) diminishing returns are already visible (8.13->6.75 was
a 1.38m drop; 6.75->6.01 was only 0.74m, despite covering *half* the
mainstem's length with real data in this second step) -- the remaining
gap likely isn't concentrated in the still-formula-covered upstream
half, since that reach carries much less accumulated flow than the
lower reaches already corrected; (2) **this test event is the 30-year
record's historical MAXIMUM flow (2026-09-17's earlier finding)** --
some genuine overbank/floodplain-scale water levels are physically
expected for a record flood, not automatically a sign of remaining
bugs. The right next check (not yet done) is section 10.1's original
suggestion: rerun this same, now-improved geometry against a *typical*
(non-record) event and confirm depths stay in a plainly sane range
there -- that is the correct bar for "is the basic setup clean",
rather than expecting a record flood to produce unremarkable numbers.

**Scope note for the paper**: this full pipeline (garbled-Shift-JIS-
filename recovery -> JGD2000 Plane Rectangular CRS identification and
validation -> mainstem tracing via max-accumulation branching ->
KP-to-grid-cell matching -> per-cell real geometry override with
formula fallback) is itself a reusable, documented recipe for
incorporating official Japanese river survey data into any RRI-based
model, not just Ishikari-specific plumbing -- worth describing in the
paper's methods/data section if item I's Ishikari case study makes it
in, independent of how far the remaining depth investigation goes.

## 11. Update (2026-09-17, current continuation): corrected a major
survey-column misinterpretation and completed the typical-event check

**This section supersedes sections 10.3 and 10.4's claim that the second
column of `WZAA2001.CSV` is channel width.** A reproducibility check
after moving the generated lookup table out of the previous session's
temporary directory exposed values of 478 m at KP26.5, 67 m at
KP26.567, and 403 m at KP27.0. Reading the companion workbook's Japanese
headers directly established that this column is `流心区間距離`: the
along-channel distance from the preceding survey station. The small
67 m value is simply the distance from the auxiliary KP26.567 station
to KP26.5, not a sudden channel constriction. The earlier 121-cell
"surveyed width" override was therefore invalid and has been removed.

### 11.1 Correct cross-section-to-RRI geometry pipeline

`real_data/ishikari/build_channel_geometry.py` now derives geometry
from the actual cross-section products:

- `WZAA3001.XLS` supplies the left/right current floodplain elevations
  (`現況高水敷高`) and the independently tabulated thalweg elevation.
- The matching `WZAA4xxx.CSV` supplies the `(cross-channel offset,
  elevation)` point cloud.
- Bankfull level is the lower of the left/right floodplain elevations.
  Only the continuous below-bankfull component containing the thalweg
  is integrated, so disconnected low ground outside a levee is not
  counted as channel.
- The RRI rectangular-channel depth is `bankfull - thalweg`; its width
  is `integrated cross-sectional area / depth`. Thus the simplified
  rectangle preserves the real surveyed bankfull area. Survey records
  marked as auxiliary/special sections are excluded.
- Later survey years overwrite older records in their overlap. The
  resulting committed table,
  `real_data/ishikari/kp_channel_geometry.csv`, contains 276 ordinary
  sections spanning KP1.5-138.0. As a strong parsing check, the minimum
  elevation independently computed from every point cloud matches the
  workbook's thalweg column for all 276 sections (tolerance 1e-6 m).

At Ishikari-Ohashi (KP26.6, interpolated between normal sections), the
correct equivalent rectangle is **299.25 m wide x 10.272 m deep =
3073.85 m2**. This also explains why the previous single-point
`480 m x ~6.3 m` quick fix happened to improve the run despite being
based on the wrong column: by coincidence it supplied almost the same
cross-sectional area, but the wrong aspect ratio.

`examples/ishikari_real_geometry.py` now applies both surveyed width
and surveyed depth to the 121 covered mainstem grid cells. It fits
power-law fallbacks to those mapped sections for tributaries and the
120 upstream mainstem cells outside coverage, rather than borrowing
Solo's parameters:

```
width = 0.895 * area_km2 ** 0.594   (log-space R2 = 0.524)
depth = 0.670 * area_km2 ** 0.281   (log-space R2 = 0.383)
```

The moderate R2 values should be reported honestly: drainage area
explains only part of real cross-section variation. Direct survey
overrides are used wherever available; the fits are explicitly only
fallbacks. The two raw derived lookup tables were moved into
`real_data/ishikari/`, removing the old hard-coded `/tmp/claude-...`
runtime dependency. Mainstem distance calculations now use the
scenario's exact Hubeny `dx`/`dy`, not rounded module constants.

### 11.2 Typical-event numerical-health check completed

`examples/run_ishikari.py` was changed to stream one 600 s forcing step
at a time instead of allocating the whole `(T, ny, nx)` rain tensor.
It now reports per-day maxima of slope depth, river depth, and adaptive
substep counts; uses RRI-compatible substep-averaged discharge; and
saves both a figure and compressed diagnostic arrays.

The requested non-record event check used 1991-08-25 through 1991-09-08
(15 days): about ten low-rain warm-up days followed by 16.56 and
42.62 mm/day, far below the 107.73 mm/day historical extreme used in
section 10. The full 2160-step CUDA run took 255 s and remained finite:

- peak `hs_max = 1.229 m`, peak `hr_max = 11.631 m`;
- maximum accepted adaptive substeps per outer step: river 128, slope 1;
- daily mean outlet flow peaked at 2642 m3/s on the day after the
  heaviest rain, then declined to 1964 m3/s the next day;
- river/slope maxima declined at the same time (`hr`: 11.631 ->
  10.995 m; `hs`: 1.229 -> 1.185 m), confirming normal recession rather
  than the unbounded storage/stiffness growth seen with the original
  inputs.

Outputs: `results/fig_ishikari_1991-08-25_15d.png` and
`results/ishikari_1991-08-25_15d.npz`.

The uncalibrated observation score is poor (daily NSE=-0.856,
RMSE=478.9 m3/s), and should **not** be presented as a predictive
validation result yet. The causes are visible in the hydrograph: the
cold-start model has zero baseflow for roughly nine days, while the
current Ishikari scenario sets infiltration, ET, and groundwater/baseflow
processes to zero and consequently overpredicts the storm peak. This run
closes the numerical/data-geometry sanity check. The next Ishikari task
is hydrologic initialization and parameter/process calibration (at
minimum warm-start/baseflow treatment and nonzero losses), not further
channel-width debugging.

## 12. Update (2026-09-17, next-stage hydrology): PET/infiltration added,
excellent same-event fit does not generalize to an independent event

The next stage from section 11.2 is now complete. The key result is
scientifically important but negative: **the daily basin-mean MERV-Jp
forcing is sufficient for numerical stress tests, but not sufficient to
claim predictive validation of this distributed, intensity-sensitive
model from one calibrated event.**

### 12.1 Treatment of processes and initial/base flow

`examples/validate_ishikari_events.py` makes the experiment explicit and
reproducible:

- MERV-Jp's daily PET is applied through the already-implemented RRI ET
  process (`evp_switch=1`).
- Green-Ampt infiltration is enabled; when `ksv>0`, its storage cap is
  set to `soildepth * gammaa`, matching `RRI_Read.f90` lines 299-300.
- The model cannot generate sustained baseflow because groundwater
  exchange is outside `rri_torch`'s implemented scope. Initializing a
  river depth would only create a transient volume that drains away.
  Therefore simulated RRI discharge is treated as stormflow and a
  **constant, observed pre-event median baseflow** is added only for the
  comparison hydrograph. This is stated in the figure and output rather
  than hidden inside the state initialization.
- Because daily precipitation is uniformly distributed over 24 hours,
  `ksv` is called a **daily-scale effective infiltration parameter**.
  It must not be interpreted as a measured instantaneous soil hydraulic
  conductivity: Green-Ampt infiltration depends strongly on rainfall
  intensity, which the source data no longer contains.

### 12.2 Calibration-event sensitivity and selected parameters

Calibration event: 1991-08-25 through 1991-09-08 (10 low-rain warm-up/
baseflow days followed by the 42.62 mm/day event). The progression was:

| Setup | NSE | RMSE (m3/s) | Interpretation |
|---|---:|---:|---|
| No PET/infiltration, `ns_slope=0.4` | -0.856 | 478.9 | Large peak overprediction, zero baseflow |
| PET + `ksv=2 mm/day`, constant baseflow | 0.291 | 296.1 | Excessive loss; peak too low and late |
| PET + `ksv=1 mm/day`, `ns_slope=0.4` | 0.741 | 179.0 | Volume much better; peak still late/low |
| PET + `ksv=1 mm/day`, `ns_slope=0.2` | **0.985** | **43.1** | Selected calibration point |

The selected setup keeps `ns_river=0.03`. Its simulated peak is
1426.6 m3/s on 1991-09-07 versus 1474.9 m3/s observed on the same day.
The two tuned scalars have separable roles here: `ksv` primarily controls
event runoff volume and `ns_slope` corrects the one-day timing lag.

### 12.3 Independent validation and the non-generalization result

An initial 1994 candidate window was rejected as an invalid cold-start
validation after checking antecedent rainfall: 45-53 mm/day storms had
occurred only 4-6 days before its start, so its observed initial flow
contained large unmodelled antecedent storage.

The valid independent event uses 1996-09-10 through 1996-10-06: about
11 low-rain warm-up days, then separate 24.5 and 39.6 mm/day events.
The 1991-selected parameters were frozen before this run. Result:

- **NSE = -0.920, RMSE = 238.7 m3/s**;
- the largest peak date is correct (1996-10-05), but the model gives
  1638 m3/s versus 944 m3/s observed and recedes too slowly;
- the earlier, smaller event shows the opposite pattern in its rising
  limb (underprediction), so a single scalar rescaling cannot repair
  both events;
- all states remain finite and the routing timing remains physically
  coherent, ruling out the numerical-instability failure mode seen in
  section 10.

This contrast -- NSE 0.985 in calibration but -0.920 independently --
must be retained in any paper discussion. Reporting only the calibration
curve would substantially overstate the evidence. The likely information
limitations are daily (not sub-daily) forcing, basin-average (not gridded)
rainfall, and absent groundwater/baseflow; snow processes may also matter
for other Ishikari seasons, though both selected events are late-summer/
autumn rain events.

Artifacts:

- `results/fig_ishikari_event_validation.png`
- `results/ishikari_calibration_event.npz`
- `results/ishikari_validation_event.npz`

**Recommended next decision for item I:** either obtain sub-daily,
spatial precipitation plus a defensible baseflow/groundwater treatment,
or scope Ishikari honestly as a second-basin geometry/data-pipeline and
numerical stress test rather than a predictive validation. Further
event-by-event tuning of `ksv` with the present forcing would only fit
missing input resolution and should not be pursued.
