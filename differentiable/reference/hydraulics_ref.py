"""Scalar, unvectorized ports of RRI_Slope.f90 (`hq`, `h2lev`) and
RRI_Riv.f90 (`hq_riv`), written to mirror the Fortran control flow
(if/elif, not torch.where) as literally as reasonable in Python.

This module has *no* dependency on rri_torch or PyTorch: it exists to be
an independent cross-check of the vectorized differentiable
implementation in `rri_torch/hydraulics.py`, not a refactor of it.
"""

from __future__ import annotations

import math


def h2lev(h: float, soildepth: float, gammaa: float) -> float:
    da_temp = soildepth * gammaa
    if soildepth == 0.0:
        return h
    elif h >= da_temp:
        return soildepth + (h - da_temp)
    else:
        rho = da_temp / soildepth
        return h / rho


def hq_slope(
    ns: float, ka: float, da: float, dm: float, beta: float,
    h: float, dh_signed: float, length: float, area: float, min_hs: float,
) -> float:
    """Mirrors RRI_Slope.f90::hq exactly, including the literal Fortran
    quirk that the `min_hs` cutoff only ever applies when `dh_signed` (as
    received at THIS call site) is positive -- see the README/module docs
    in rri_torch for why the differentiable model instead applies it
    symmetrically. With `min_wc4latflow = 0` (as in every scenario this
    reference is exercised against) `min_hs` is always 0 and this
    distinction never actually fires.
    """
    if dh_signed > 0.0 and h <= min_hs:
        q = 0.0
    else:
        dh = abs(dh_signed)
        km = ka / beta if beta > 0.0 else 0.0
        vm = km * dh
        va = ka * dh if da > 0.0 else 0.0
        al = math.sqrt(dh) / ns
        m = 5.0 / 3.0

        if h < dm:
            q = vm * dm * (h / dm) ** beta
        elif h < da:
            q = vm * dm + va * (h - dm)
        else:
            q = vm * dm + va * (h - dm) + al * (h - da) ** m

    return q * length / area


def hq_river(h: float, dh_signed: float, width: float, ns_river: float) -> float:
    dh = abs(dh_signed)
    a = math.sqrt(dh) / ns_river
    m = 5.0 / 3.0
    return a * (h ** m) * width
