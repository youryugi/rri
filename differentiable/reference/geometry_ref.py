"""Independent (no rri_torch import) grid/river-tree construction, mirroring
RRI_Sub.f90's index building but operating directly on (i, j) grids with
plain Python loops instead of compacted 1-D index arrays.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np

D8_OFFSETS = {
    1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1),
    16: (0, -1), 32: (-1, -1), 64: (-1, 0), 128: (-1, 1),
}

# 8-direction diffusive-wave slope edges, matching RRI_Slope.f90's
# l = 1..4 loop (E, S, SE, SW) and RRI_Sub.f90's l1=dy/2, l2=dx/2, l3=diag/4.
SLOPE_EDGES = ["E", "S", "SE", "SW"]
SLOPE_EDGE_OFFSETS = {"E": (0, 1), "S": (1, 0), "SE": (1, 1), "SW": (1, -1)}


@dataclass
class RefGrid:
    ny: int
    nx: int
    dx: float
    dy: float
    area: float = field(init=False)
    domain: np.ndarray = None
    sink: np.ndarray = None
    riv_mask: np.ndarray = None
    zb: np.ndarray = None
    down_riv: dict = None      # (i, j) -> (ii, jj) or None
    dis_riv: dict = None       # (i, j) -> distance [m]
    len_riv_grid: np.ndarray = None

    def __post_init__(self):
        self.area = self.dx * self.dy

    def slope_edge_geometry(self, edge: str):
        diag = math.hypot(self.dx, self.dy)
        return {
            "E": (self.dx, self.dy / 2.0),
            "S": (self.dy, self.dx / 2.0),
            "SE": (diag, diag / 4.0),
            "SW": (diag, diag / 4.0),
        }[edge]

    def neighbors_slope(self, i, j):
        """Yields (edge, ii, jj, distance, length) for the 4 outgoing
        half-edges of cell (i, j) that stay inside the active domain."""
        for edge in SLOPE_EDGES:
            di, dj = SLOPE_EDGE_OFFSETS[edge]
            ii, jj = i + di, j + dj
            if not (0 <= ii < self.ny and 0 <= jj < self.nx):
                continue
            if not (self.domain[i, j] and self.domain[ii, jj]):
                continue
            distance, length = self.slope_edge_geometry(edge)
            yield edge, ii, jj, distance, length

    def inflow_neighbors_slope(self, i, j):
        """Yields (edge, ii, jj) of every neighbour whose outgoing `edge`
        points INTO (i, j) -- the mirror image of `neighbors_slope`,
        replicating what `up_slo_gather` collects."""
        for edge in SLOPE_EDGES:
            di, dj = SLOPE_EDGE_OFFSETS[edge]
            ii, jj = i - di, j - dj
            if not (0 <= ii < self.ny and 0 <= jj < self.nx):
                continue
            if not (self.domain[i, j] and self.domain[ii, jj]):
                continue
            yield edge, ii, jj


def build_ref_grid(
    domain: np.ndarray,
    zb: np.ndarray,
    dir_grid: np.ndarray,
    riv_mask: np.ndarray,
    len_riv: np.ndarray,
    dx: float,
    dy: float,
    sink: np.ndarray | None = None,
) -> RefGrid:
    ny, nx = domain.shape
    domain_b = domain.astype(bool)
    riv_b = riv_mask.astype(bool) & domain_b

    down_riv = {}
    dis_riv = {}
    auto_sink = np.zeros((ny, nx), dtype=bool)
    for i in range(ny):
        for j in range(nx):
            if not riv_b[i, j]:
                continue
            d = int(dir_grid[i, j])
            if d in D8_OFFSETS:
                di, dj = D8_OFFSETS[d]
                ii, jj = i + di, j + dj
                if 0 <= ii < ny and 0 <= jj < nx and riv_b[ii, jj]:
                    down_riv[(i, j)] = (ii, jj)
                    dist = math.hypot(dx, dy) if (di != 0 and dj != 0) else (dy if di != 0 else dx)
                    dis_riv[(i, j)] = 0.5 * (len_riv[i, j] + len_riv[ii, jj])
                    continue
            down_riv[(i, j)] = None
            auto_sink[i, j] = True
            dist = math.hypot(dx, dy) if (int(dir_grid[i, j]) in D8_OFFSETS and
                                           D8_OFFSETS[int(dir_grid[i, j])][0] != 0 and
                                           D8_OFFSETS[int(dir_grid[i, j])][1] != 0) else dx
            dis_riv[(i, j)] = 0.5 * (len_riv[i, j] + dist)

    sink_b = auto_sink if sink is None else sink.astype(bool)

    return RefGrid(
        ny=ny, nx=nx, dx=dx, dy=dy,
        domain=domain_b, sink=sink_b, riv_mask=riv_b, zb=zb.astype(np.float64),
        down_riv=down_riv, dis_riv=dis_riv, len_riv_grid=len_riv.astype(np.float64),
    )
