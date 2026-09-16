"""Static grid geometry: domain/river masks, elevations, and the river tree.

This mirrors the index-building code in RRI_Sub.f90 (`sub_slo_ij2idx`,
`sub_riv_ij2idx`, and the D8 downstream search) but is built once with plain
Python/NumPy -- none of it needs to be differentiable, only the physics that
runs on top of it does.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np
import torch

# RRI / ESRI D8 flow-direction codes -> (di, dj) neighbour offset.
# 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE (see RRI_Sub.f90).
D8_OFFSETS = {
    1: (0, 1),
    2: (1, 1),
    4: (1, 0),
    8: (1, -1),
    16: (0, -1),
    32: (-1, -1),
    64: (-1, 0),
    128: (-1, 1),
}


@dataclass
class Grid:
    """Static (non-learnable) description of a catchment.

    All tensors use shape (ny, nx) for grid fields and (n_riv,) for
    river-cell fields, matching the *_idx 1-D compaction RRI performs
    internally for river cells.
    """

    ny: int
    nx: int
    dx: float
    dy: float
    domain: torch.Tensor      # bool (ny, nx): cell participates in the model
    sink: torch.Tensor        # bool (ny, nx): outlet cell, drains to zero each step (domain==2 in RRI)
    riv_mask: torch.Tensor    # bool (ny, nx): cell contains a river channel
    zb: torch.Tensor          # float (ny, nx): ground (slope) elevation [m]

    # river-cell 1-D arrays
    riv_i: torch.Tensor       # long (n_riv,)
    riv_j: torch.Tensor       # long (n_riv,)
    riv_index_map: torch.Tensor  # long (ny, nx): -1 if not a river cell
    down_riv: torch.Tensor    # long (n_riv,): index of downstream river cell, -1 if outlet
    len_riv: torch.Tensor     # float (n_riv,): channel reach length in this cell [m]
    dis_riv: torch.Tensor     # float (n_riv,): distance to downstream cell [m]

    area: float = field(init=False)

    def __post_init__(self):
        self.area = self.dx * self.dy

    @property
    def n_riv(self) -> int:
        return self.riv_i.shape[0]

    def to(self, *args, **kwargs) -> "Grid":
        tensor_fields = [
            "domain", "sink", "riv_mask", "zb", "riv_i", "riv_j",
            "riv_index_map", "down_riv", "len_riv", "dis_riv",
        ]
        kwargs2 = dict(kwargs)
        for f in tensor_fields:
            setattr(self, f, getattr(self, f).to(*args, **kwargs2))
        return self

    @classmethod
    def build(
        cls,
        domain: np.ndarray,
        zb: np.ndarray,
        dir_grid: np.ndarray,
        riv_mask: np.ndarray | None = None,
        sink: np.ndarray | None = None,
        len_riv: np.ndarray | None = None,
        dx: float = 1.0,
        dy: float = 1.0,
        dtype: torch.dtype = torch.float64,
    ) -> "Grid":
        """Build a Grid from plain NumPy arrays.

        Parameters
        ----------
        domain : bool (ny, nx), True where the model is active.
        zb : float (ny, nx), ground elevation.
        dir_grid : int (ny, nx), D8 flow direction code (RRI convention),
            used only to trace the river network's downstream tree.
        riv_mask : bool (ny, nx), True where a river channel exists.
            Defaults to all-False (pure overland-flow model).
        sink : bool (ny, nx), outlet cells that drain to zero every step.
            Defaults to river cells with no valid downstream river cell.
        len_riv : float (ny, nx) channel length within each river cell.
            Defaults to dx for every river cell.
        """
        ny, nx = domain.shape
        domain_b = domain.astype(bool)
        if riv_mask is None:
            riv_mask_b = np.zeros_like(domain_b)
        else:
            riv_mask_b = riv_mask.astype(bool) & domain_b

        riv_cells = list(zip(*np.where(riv_mask_b)))
        n_riv = len(riv_cells)
        riv_index_map = -np.ones((ny, nx), dtype=np.int64)
        for n, (i, j) in enumerate(riv_cells):
            riv_index_map[i, j] = n

        if len_riv is None:
            len_riv_grid = np.full((ny, nx), dx, dtype=np.float64)
        else:
            len_riv_grid = len_riv

        down_riv = -np.ones(n_riv, dtype=np.int64)
        dis_riv = np.zeros(n_riv, dtype=np.float64)
        riv_i = np.zeros(n_riv, dtype=np.int64)
        riv_j = np.zeros(n_riv, dtype=np.int64)
        auto_sink = np.zeros((ny, nx), dtype=bool)
        for n, (i, j) in enumerate(riv_cells):
            riv_i[n] = i
            riv_j[n] = j
            d = int(dir_grid[i, j])
            down_len = len_riv_grid[i, j]
            if d in D8_OFFSETS:
                di, dj = D8_OFFSETS[d]
                ii, jj = i + di, j + dj
                if 0 <= ii < ny and 0 <= jj < nx and riv_index_map[ii, jj] >= 0:
                    down_riv[n] = riv_index_map[ii, jj]
                    down_len = len_riv_grid[ii, jj]
                    dist = math.hypot(dx, dy) if (di != 0 and dj != 0) else (dy if di != 0 else dx)
                else:
                    auto_sink[i, j] = True
                    dist = math.hypot(dx, dy) if (di != 0 and dj != 0) else (dy if di != 0 else dx)
            else:
                auto_sink[i, j] = True
                dist = dx
            dis_riv[n] = 0.5 * (len_riv_grid[i, j] + down_len) if down_riv[n] >= 0 else 0.5 * (len_riv_grid[i, j] + dist)

        if sink is None:
            sink_b = auto_sink
        else:
            sink_b = sink.astype(bool)

        def T(a, kind=dtype):
            return torch.as_tensor(a, dtype=kind)

        return cls(
            ny=ny, nx=nx, dx=dx, dy=dy,
            domain=T(domain_b, torch.bool),
            sink=T(sink_b, torch.bool),
            riv_mask=T(riv_mask_b, torch.bool),
            zb=T(zb),
            riv_i=T(riv_i, torch.long),
            riv_j=T(riv_j, torch.long),
            riv_index_map=T(riv_index_map, torch.long),
            down_riv=T(down_riv, torch.long),
            len_riv=T(len_riv_grid[riv_i, riv_j] if n_riv else np.zeros(0)),
            dis_riv=T(dis_riv),
        )

    def gather_to_riv(self, grid_field: torch.Tensor) -> torch.Tensor:
        """(ny, nx) -> (n_riv,) at river-cell locations."""
        return grid_field[self.riv_i, self.riv_j]

    def scatter_from_riv(self, riv_field: torch.Tensor, fill_value: float = 0.0) -> torch.Tensor:
        """(n_riv,) -> (ny, nx), zero (or fill_value) off the river cells."""
        out = torch.full((self.ny, self.nx), fill_value, dtype=riv_field.dtype, device=riv_field.device)
        if self.n_riv:
            out[self.riv_i, self.riv_j] = riv_field
        return out
