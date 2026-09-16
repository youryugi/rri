"""First real-data validation run: the Solo River basin (see
differentiable/HANDOFF.md section 3-5 for full background).

Runs `rri_torch` on the real 336x204 / 18582-active-cell Solo catchment
and reports outlet discharge at all 6 gauge stations, plus basic
stability/mass-balance diagnostics. Meant to be run incrementally:
start with a short time window (`--steps`) to catch stability/NaN issues
cheaply before committing to the full 360h/2160-step event.

Usage:
  python examples/solo_validation.py --steps 20 --substeps 20
  python examples/solo_validation.py --steps 2160 --substeps 20 --device cuda
"""
from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch

from rri_torch.model import RRIModel
from real_data.solo_river import build_solo_scenario, build_rain_sequence, read_esri_ascii

PROJECT_DIR = "/home/yang/github/rri/RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s"
DT = 600.0  # seconds, matches RRI_Input.txt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=20, help="number of dt=600s outer steps to run")
    ap.add_argument("--substeps-slope", type=int, default=20)
    ap.add_argument("--substeps-river", type=int, default=20)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--save", default=None, help="path to save qr_outlet series + gauge names as .pt")
    args = ap.parse_args()

    device = torch.device(args.device)
    dtype = torch.float64

    print(f"[load] building Solo scenario from {PROJECT_DIR} ...")
    t0 = time.time()
    grid, params, outlets, rain = build_solo_scenario(PROJECT_DIR, dtype=dtype)
    print(f"[load] done in {time.time()-t0:.1f}s: ny={grid.ny} nx={grid.nx} "
          f"n_riv={grid.n_riv} active={int(grid.domain.sum())} dx={grid.dx:.2f} dy={grid.dy:.2f}")
    print(f"[load] outlets: {outlets}")

    dem = read_esri_ascii(os.path.join(PROJECT_DIR, "topo", "adem.txt"))
    print(f"[load] building rain sequence for {args.steps} steps ...")
    t0 = time.time()
    rain_seq = build_rain_sequence(
        rain, grid.nx, grid.ny, dem.xllcorner, dem.yllcorner, dem.cellsize,
        DT, args.steps, dtype=dtype,
    )
    print(f"[load] rain sequence built in {time.time()-t0:.1f}s, "
          f"peak rate = {rain_seq.max().item()*3600*1000:.3f} mm/h")

    grid = grid.to(device)
    params.slope.ns = params.slope.ns.to(device)
    params.slope.ka = params.slope.ka.to(device)
    params.slope.beta = params.slope.beta.to(device)
    params.slope.da = params.slope.da.to(device)
    params.slope.dm = params.slope.dm.to(device)
    params.slope.soildepth = params.slope.soildepth.to(device)
    params.slope.gammaa = params.slope.gammaa.to(device)
    params.slope.min_wc4latflow = params.slope.min_wc4latflow.to(device)
    params.river.width = params.river.width.to(device)
    params.river.depth = params.river.depth.to(device)
    params.river.ns_river = params.river.ns_river.to(device)
    params.river.height = params.river.height.to(device)
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)
    rain_seq = rain_seq.to(device)

    model = RRIModel(grid, n_substeps_slope=args.substeps_slope, n_substeps_river=args.substeps_river)

    # `diag.qr[outlet_riv_index]` supports fancy indexing with a list, so one
    # simulate() call reports every gauge's hydrograph.
    outlet_names = list(outlets.keys())
    outlet_idx = [outlets[n] for n in outlet_names]

    init_state = tuple(x.to(device) for x in model.initial_state(dtype=dtype))

    print(f"[run] simulating {args.steps} steps of dt={DT}s on {device} ...")
    t0 = time.time()
    out = model.simulate(rain_seq, params, DT, outlet_riv_index=outlet_idx,
                          init_state=init_state, record_full_state=False)
    if device.type == "cuda":
        torch.cuda.synchronize()
    print(f"[run] done in {time.time()-t0:.1f}s")

    hs_final, hr_final = out["hs_final"], out["hr_final"]
    print(f"[check] hs finite: {torch.isfinite(hs_final).all().item()}, "
          f"max hs = {hs_final.max().item():.4f} m")
    print(f"[check] hr finite: {torch.isfinite(hr_final).all().item()}, "
          f"max hr = {hr_final.max().item():.4f} m")

    qr_all = out["qr_outlet"]  # (T, n_gauges)
    for k, name in enumerate(outlet_names):
        qr = qr_all[:, k]
        print(f"[gauge {name:10s}] peak qr = {qr.max().item():10.3f} m^3/s "
              f"at step {int(qr.argmax())}/{args.steps}, finite={torch.isfinite(qr).all().item()}")

    if args.save:
        torch.save({"names": outlet_names, "qr": qr_all.cpu(), "dt": DT, "steps": args.steps}, args.save)
        print(f"[save] wrote {args.save}")


if __name__ == "__main__":
    main()
