"""Small vectorized helpers shared by the slope and river routing modules.

The RRI slope scheme is a structured-grid, 8-direction diffusive wave. Rather
than reproducing the Fortran's explicit neighbour-index arrays
(``down_slo_idx`` / ``up_slo_gather_*`` in ``RRI_Slope.f90``), we exploit the
fact that our state lives on a regular (ny, nx) grid and get every neighbour
via a single padded shift. This keeps the whole slope update as a handful of
tensor ops, which is what makes it cheap to differentiate through.
"""

from __future__ import annotations

import torch
import torch.nn.functional as F

# The four "forward" directions used by RRI's 8-direction diffusive scheme.
# Every one of the 8 neighbours of a cell is reached either directly through
# one of these four offsets, or as the mirror image of the same offset
# applied at the neighbouring cell (see `neighbor_of` below for the inflow
# side). This matches `down_slo_idx(l, k)` for l = 1..4 in RRI_Slope.f90:
#   l=1 -> East, l=2 -> South, l=3 -> South-East, l=4 -> South-West
EDGE_OFFSETS = {
    "E": (0, 1),
    "S": (1, 0),
    "SE": (1, 1),
    "SW": (1, -1),
}


def shift(t: torch.Tensor, di: int, dj: int, pad_value: float = 0.0) -> torch.Tensor:
    """Return an array `out` with out[..., i, j] = t[..., i+di, j+dj].

    Cells that would fall outside the grid are filled with `pad_value`.
    `di`, `dj` must be in {-1, 0, 1} (all we need for 8-connected grids).
    """
    assert -1 <= di <= 1 and -1 <= dj <= 1
    ny, nx = t.shape[-2], t.shape[-1]
    padded = F.pad(t, (1, 1, 1, 1), value=pad_value)
    return padded[..., 1 + di:1 + di + ny, 1 + dj:1 + dj + nx]


def safe_sqrt(x: torch.Tensor, eps: float = 1e-12) -> torch.Tensor:
    """sqrt(x) with the argument floored at `eps` before the square root.

    Plain `torch.sqrt(torch.clamp(x, min=0.0))` is forward-safe but not
    backward-safe: d(sqrt)/dx = 1/(2*sqrt(x)) diverges as x -> 0, and even
    when that branch is masked out by `torch.where` at some downstream
    point, autograd still multiplies the (infinite) local gradient by a
    zero mask -- 0 * inf = NaN. Flooring at a small positive `eps` instead
    of exactly 0 keeps the local gradient finite everywhere, at the cost
    of a negligible (~sqrt(eps)) forward-value bias for arguments that are
    exactly zero.
    """
    return torch.sqrt(torch.clamp(x, min=eps))


def inflow_shift(t: torch.Tensor, edge: str, pad_value: float = 0.0) -> torch.Tensor:
    """For a per-cell outgoing-edge quantity `t` (e.g. flux leaving each cell
    toward `edge`), return the value contributed *into* each cell by its
    neighbour on the opposite side of that same edge.

    e.g. inflow_shift(qE, "E") gives, at (i, j), the East-directed flux
    leaving the cell to the West, i.e. qE[i, j-1] -- exactly the term
    `up_slo_gather` collects in the Fortran implementation.
    """
    di, dj = EDGE_OFFSETS[edge]
    # the neighbour that flows INTO (i, j) via this edge sits at (i-di, j-dj)
    return shift(t, -di, -dj, pad_value=pad_value)
