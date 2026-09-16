"""Diagnoses the Cepu-gauge backwater-like divergence found in
compare_solo_hydro.py (see HANDOFF.md section 6a): runs the full 360h
Solo event while recording discharge at every river cell along the chain
from Cepu down to the true basin outlet (147 cells), so we can see
whether the anomalous growth originates right at the outlet's free-flow
boundary and propagates upstream, or somewhere else in the reach.

Also supports overriding `outlet_slope` (default in RiverParams is 1e-3;
the real local bed slope over the last 5 outlet-reach cells was measured
at ~2.2e-4 -- see HANDOFF.md) to test whether that mismatch is the cause.

Usage:
  python examples/solo_outlet_diagnosis.py --outlet-slope 1e-3 --save diag_default.pt
  python examples/solo_outlet_diagnosis.py --outlet-slope 2.2e-4 --save diag_fixedslope.pt
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
DT = 600.0


def chain_from(grid, start_idx):
    idx = start_idx
    chain = [idx]
    while grid.down_riv[idx].item() >= 0:
        idx = grid.down_riv[idx].item()
        chain.append(idx)
    return chain


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=2160)
    ap.add_argument("--substeps-slope", type=int, default=40)
    ap.add_argument("--substeps-river", type=int, default=40)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--outlet-slope", type=float, default=None,
                     help="override RiverParams.outlet_slope (default in the model is 1e-3)")
    ap.add_argument("--save", required=True)
    args = ap.parse_args()

    device = torch.device(args.device)
    dtype = torch.float64

    grid, params, outlets, rain = build_solo_scenario(PROJECT_DIR, dtype=dtype)
    cepu = outlets["Cepu"]
    chain = chain_from(grid, cepu)  # Cepu ... true outlet, 148 entries
    print(f"[chain] Cepu -> outlet: {len(chain)} cells")

    if args.outlet_slope is not None:
        params.river.outlet_slope = args.outlet_slope
        print(f"[param] outlet_slope overridden to {args.outlet_slope}")
    else:
        print(f"[param] outlet_slope left at default {params.river.outlet_slope}")

    dem = read_esri_ascii(os.path.join(PROJECT_DIR, "topo", "adem.txt"))
    rain_seq = build_rain_sequence(
        rain, grid.nx, grid.ny, dem.xllcorner, dem.yllcorner, dem.cellsize,
        DT, args.steps, dtype=dtype,
    )

    grid = grid.to(device)
    for obj, fields in [
        (params.slope, ["ns", "ka", "beta", "da", "dm", "soildepth", "gammaa", "min_wc4latflow"]),
        (params.river, ["width", "depth", "ns_river", "height"]),
    ]:
        for f in fields:
            setattr(obj, f, getattr(obj, f).to(device))
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)
    rain_seq = rain_seq.to(device)

    model = RRIModel(grid, n_substeps_slope=args.substeps_slope, n_substeps_river=args.substeps_river)
    init_state = tuple(x.to(device) for x in model.initial_state(dtype=dtype))

    # Record qr along the whole Cepu->outlet chain (148 points) -- cheap,
    # it's just fancy-indexing the existing per-step qr diagnostic.
    print(f"[run] simulating {args.steps} steps on {device} ...")
    t0 = time.time()
    out = model.simulate(rain_seq, params, DT, outlet_riv_index=chain,
                          init_state=init_state, record_full_state=False)
    if device.type == "cuda":
        torch.cuda.synchronize()
    print(f"[run] done in {time.time()-t0:.1f}s")

    qr_chain = out["qr_outlet"].cpu()  # (T, len(chain))
    finite = torch.isfinite(qr_chain).all(dim=0)
    print(f"[check] finite along chain: {finite.all().item()} "
          f"({(~finite).sum().item()} non-finite cells)")

    torch.save({
        "chain": chain,
        "chain_i": grid.riv_i.cpu()[chain].tolist(),
        "chain_j": grid.riv_j.cpu()[chain].tolist(),
        "qr_chain": qr_chain,
        "dt": DT,
        "steps": args.steps,
        "outlet_slope": params.river.outlet_slope,
    }, args.save)
    print(f"[save] wrote {args.save}")

    # Quick summary: peak qr and peak step for a handful of positions
    # along the chain (0=Cepu, -1=true outlet).
    positions = [0, len(chain) // 4, len(chain) // 2, 3 * len(chain) // 4,
                 len(chain) - 15, len(chain) - 5, len(chain) - 1]
    for p in sorted(set(positions)):
        qr = qr_chain[:, p]
        print(f"[chain pos {p:4d}/{len(chain)-1}] peak qr = {qr.max().item():10.3f} m^3/s "
              f"at step {int(qr.argmax())}/{args.steps}")


if __name__ == "__main__":
    main()
