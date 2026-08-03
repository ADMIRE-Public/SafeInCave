# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Plain-NumPy Mandel-vector post-processing helpers shared by dolfinx-based
output classes — distinct from any JAX-traced, differentiable math consumed
directly inside a return-mapping solve.
"""

from __future__ import annotations

import numpy as np

SQRT2 = np.sqrt(2.0)


def mandel_to_tensor3x3(v: np.ndarray) -> np.ndarray:
    """Batched Mandel vector(s) (N, 4|6) -> full tensors (N, 3, 3).

    4-component vectors are plane strain: [xx, yy, zz, sqrt2*xy].
    """
    v = np.atleast_2d(v)
    n, dim = v.shape
    T = np.zeros((n, 3, 3))
    T[:, 0, 0] = v[:, 0]
    T[:, 1, 1] = v[:, 1]
    T[:, 2, 2] = v[:, 2]
    T[:, 0, 1] = T[:, 1, 0] = v[:, 3] / SQRT2
    if dim == 6:
        T[:, 0, 2] = T[:, 2, 0] = v[:, 4] / SQRT2
        T[:, 1, 2] = T[:, 2, 1] = v[:, 5] / SQRT2
    return T


def p_q_from_mandel(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean stress (tr/3, signed) and von Mises stress from Mandel vectors."""
    v = np.atleast_2d(v)
    sm = (v[:, 0] + v[:, 1] + v[:, 2]) / 3.0
    dev = v.copy()
    dev[:, :3] -= sm[:, None]
    # Mandel components carry the sqrt(2) factor: |dev|^2 is the tensor norm
    q = np.sqrt(1.5 * np.sum(dev**2, axis=1))
    return sm, q


def principal_from_tensor(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Batched eigen-decomposition; eigenvalues descending (sig1 >= sig3)."""
    w, V = np.linalg.eigh(T)
    return w[:, ::-1], V[:, :, ::-1]
