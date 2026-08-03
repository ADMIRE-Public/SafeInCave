# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Quadrature-space view helper shared by dolfinx-based output classes
(e.g. the autodiff/JAX external-operator solver path's ``SaveFields``
and ``PointLogger``).
"""

from __future__ import annotations

import numpy as np
from dolfinx import fem


class _QuadratureView:
    """Cell-wise view of a quadrature-space Function via its dofmap."""

    def __init__(self, fn: fem.Function, n_comp: int):
        self.fn = fn
        self.n_comp = n_comp
        self.dm = fn.function_space.dofmap.list  # (n_cells, n_qp)
        self.n_cells, self.n_qp = self.dm.shape

    def cell_values(self) -> np.ndarray:
        """(n_cells, n_qp, n_comp) array of the current values."""
        flat = self.fn.x.array.reshape(-1, self.n_comp)
        return flat[self.dm]

    def cell_mean(self) -> np.ndarray:
        """(n_cells, n_comp) quadrature-point average per cell."""
        return self.cell_values().mean(axis=1)
