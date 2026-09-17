"""Direct comparison against the real compiled Fortran binary's own
spatial-field output (`out/qr_*.out`), not just the single-point
`hydro.txt` -- see HANDOFF.md section 6a for how the native build was
finally gotten working (an actual bug in RRI.f90's rain.dat-reading loop,
not a performance quirk).

Usage:
  python examples/compare_solo_spatial.py --run <solo_validation.py --save output>
"""
from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import torch

from real_data.solo_river import read_locations

PROJECT_DIR = "/home/yang/github/rri/RRI_1_4_2_7_GUI_Beta/RRI-CUI/Project/solo30s"
FORTRAN_OUT_DIR = "/tmp/claude-1004/-home-yang-github-rri/56b4d09e-3c12-4eb2-94cb-037713c8deea/scratchpad/solo30s_run/out"
MAXT, OUTNUM = 2160, 96


def fortran_nint(x: float) -> int:
    return math.floor(x + 0.5) if x >= 0 else math.ceil(x - 0.5)


def output_step_indices() -> list[int]:
    """The 1-indexed outer-step `t` at which each of the 96 `out/*_NNNNNN.out`
    files was written, exactly matching RRI.f90's `out_next = nint((tt+1) *
    dble(maxt)/dble(outnum))` (lines 593-1013)."""
    out_dt = MAXT / OUTNUM
    return [fortran_nint(k * out_dt) for k in range(1, OUTNUM + 1)]


def read_grid_out(path: str, ny: int, nx: int) -> np.ndarray:
    data = np.loadtxt(path)
    assert data.shape == (ny, nx), f"{path}: expected ({ny},{nx}), got {data.shape}"
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="a solo_outlet_diagnosis.py --save .pt file "
                                                   "(needs qr_chain covering the cells of interest)")
    ap.add_argument("--instantaneous", action="store_true",
                     help="use qr_chain (instantaneous) even if qr_avg_chain was saved")
    args = ap.parse_args()

    d = torch.load(args.run, weights_only=False)
    chain = d["chain"]
    chain_i, chain_j = d["chain_i"], d["chain_j"]
    use_avg = ("qr_avg_chain" in d) and not args.instantaneous
    qr_chain = (d["qr_avg_chain"] if use_avg else d["qr_chain"]).numpy()  # (T, len(chain))
    print(f"[sim] comparing against {'qr_avg (time-averaged)' if use_avg else 'qr (instantaneous)'}")
    T = qr_chain.shape[0]
    print(f"[sim] loaded {T}-step run, chain of {len(chain)} cells")

    ny, nx = 204, 336
    steps = output_step_indices()
    print(f"[fortran] output files map to outer steps: first={steps[0]} last={steps[-1]}")

    # Compare at every 8th chain position (0=Cepu ... last=outlet) across
    # every 8th output file, to keep this readable while still covering
    # the whole reach and the whole event. The true outlet cell itself
    # (last position) is excluded from stats: RRI_Riv.f90's qr_calc skips
    # any cell with domain_riv_idx==2 entirely ("if(domain_riv_idx(k).eq.2)
    # cycle") -- the sink cell's own qr just stays at its 0 initialization
    # in the real output, whereas rri_torch's river_rhs does compute an
    # actual free-outflow discharge there. Different reporting convention
    # at that one cell, not a physics disagreement -- comparing it against
    # a near-zero reference blows up into meaningless percentages.
    outlet_pos = len(chain) - 1
    positions = sorted(set(range(0, len(chain), 8)) | {outlet_pos})
    file_ks = list(range(8, OUTNUM + 1, 8))

    header = "step/pos".ljust(10) + "".join(f"pos{p:4d}".rjust(11) for p in positions)
    print(f"\n=== qr: rri_torch vs real Fortran out/qr_*.out (pos {outlet_pos} = outlet, excluded from stats) ===")
    print(header)
    all_reldiffs = []
    for k in file_ks:
        t_step = steps[k - 1]
        if t_step > T:
            continue
        fpath = os.path.join(FORTRAN_OUT_DIR, f"qr_{k:06d}.out")
        grid = read_grid_out(fpath, ny, nx)
        row = f"t={t_step:5d}".ljust(10)
        for p in positions:
            i, j = chain_i[p], chain_j[p]
            f_val = grid[i, j]
            t_val = qr_chain[t_step - 1, p]
            if f_val < 0:  # -0.1 sentinel = outside domain / not a river cell in fortran's own mask
                row += "     n/a   "
                continue
            reldiff = (t_val - f_val) / max(abs(f_val), 1.0) * 100
            row += f"{reldiff:+9.1f}% "
            if p != outlet_pos:
                all_reldiffs.append(reldiff)
        print(row)

    if all_reldiffs:
        arr = np.array(all_reldiffs)
        print(f"\n[summary] %diff over {len(arr)} (step,cell) samples, outlet cell excluded: "
              f"mean={arr.mean():+.1f}%  median={np.median(arr):+.1f}%  "
              f"|mean|={np.abs(arr).mean():.1f}%  std={arr.std():.1f}%")


if __name__ == "__main__":
    main()
