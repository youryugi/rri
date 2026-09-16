"""2D diffusive-wave slope (overland + subsurface) routing.

Vectorized re-implementation of `funcs` / `qs_calc_do` in RRI_Slope.f90,
restricted to the "8-direction" diffusive-wave mode (`eight_dir = 1`),
which is RRI's standard/default slope scheme. The single-direction
kinematic-wave mode (`eight_dir = 0`, `dif = 0` land use) is not
implemented -- see README for scope notes.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import torch

from .hydraulics import h2lev, hq_slope
from .tensor_ops import shift, inflow_shift, EDGE_OFFSETS

# Edge -> (distance factor, cross-flow length factor), as fractions of
# (dx, dy). Matches RRI_Sub.f90: l1 = dy/2, l2 = dx/2, l3 = sqrt(dx^2+dy^2)/4
# for the 8-direction scheme.
_EDGES = ("E", "S", "SE", "SW")


@dataclass
class SlopeParams:
    """Per-cell (ny, nx) slope parameters, all learnable in principle."""

    ns: torch.Tensor              # Manning's roughness [s/m^(1/3)]
    ka: torch.Tensor              # lateral saturated hydraulic conductivity [m/s]
    beta: torch.Tensor            # unsaturated-flow exponent [-]
    da: torch.Tensor              # saturated-storage depth threshold [m]
    dm: torch.Tensor              # unsaturated/linear transition depth [m]
    soildepth: torch.Tensor       # soil layer thickness [m]
    gammaa: torch.Tensor          # effective porosity [-]
    min_wc4latflow: torch.Tensor  # min water content fraction before lateral flow starts [-]


def _edge_geometry(dx: float, dy: float):
    diag = math.hypot(dx, dy)
    return {
        "E": (dx, dy / 2.0),
        "S": (dy, dx / 2.0),
        "SE": (diag, diag / 4.0),
        "SW": (diag, diag / 4.0),
    }


def slope_edge_fluxes(
    hs: torch.Tensor,
    zb: torch.Tensor,
    domain: torch.Tensor,
    params: SlopeParams,
    dx: float,
    dy: float,
    area: float,
) -> dict[str, torch.Tensor]:
    """Return {edge_name: outflow flux [m/s], positive = leaving the cell}
    for each of the 4 half-edges (E, S, SE, SW), following RRI's
    diffusive-wave sign convention.
    """
    geom = _edge_geometry(dx, dy)
    lev = h2lev(hs, params.soildepth, params.gammaa)
    fluxes = {}

    for edge in _EDGES:
        di, dj = EDGE_OFFSETS[edge]
        distance, length = geom[edge]

        zb_n = shift(zb, di, dj)
        hs_n = shift(hs, di, dj)
        lev_n = shift(lev, di, dj)
        domain_n = shift(domain, di, dj, pad_value=False)
        valid = domain & domain_n

        ns_n = shift(params.ns, di, dj)
        ka_n = shift(params.ka, di, dj)
        beta_n = shift(params.beta, di, dj)
        da_n = shift(params.da, di, dj)
        dm_n = shift(params.dm, di, dj)
        soildepth_n = shift(params.soildepth, di, dj)
        gammaa_n = shift(params.gammaa, di, dj)
        min_wc_n = shift(params.min_wc4latflow, di, dj)

        dh = ((zb + lev) - (zb_n + lev_n)) / distance

        hw_out = torch.where(zb >= zb_n, hs, torch.clamp(zb + hs - zb_n, min=0.0))
        hw_in = torch.where(zb_n >= zb, hs_n, torch.clamp(zb_n + hs_n - zb, min=0.0))

        q_out = hq_slope(
            params.ns, params.ka, params.da, params.dm, params.beta,
            hw_out, dh, length, area,
            params.soildepth, params.gammaa, params.min_wc4latflow,
        )
        q_in = hq_slope(
            ns_n, ka_n, da_n, dm_n, beta_n,
            hw_in, dh, length, area,
            soildepth_n, gammaa_n, min_wc_n,
        )

        q_edge = torch.where(dh >= 0, q_out, -q_in)
        fluxes[edge] = q_edge * valid
    return fluxes


def slope_rhs(
    hs: torch.Tensor,
    zb: torch.Tensor,
    domain: torch.Tensor,
    params: SlopeParams,
    rain: torch.Tensor,
    dx: float,
    dy: float,
    area: float,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    """d(hs)/dt [m/s] and the per-edge fluxes (for diagnostics / the
    river-slope exchange step), following `funcs` in RRI_Slope.f90.

    Note: unlike the Fortran, infiltration and evapotranspiration are NOT
    included here either -- RRI applies them as a separate operator-split
    update after routing completes (see `model.py`), and we preserve that
    splitting.
    """
    fluxes = slope_edge_fluxes(hs, zb, domain, params, dx, dy, area)
    outflow = fluxes["E"] + fluxes["S"] + fluxes["SE"] + fluxes["SW"]
    inflow = (
        inflow_shift(fluxes["E"], "E")
        + inflow_shift(fluxes["S"], "S")
        + inflow_shift(fluxes["SE"], "SE")
        + inflow_shift(fluxes["SW"], "SW")
    )
    dhsdt = rain - outflow + inflow
    dhsdt = torch.where(domain, dhsdt, torch.zeros_like(dhsdt))
    return dhsdt, fluxes
