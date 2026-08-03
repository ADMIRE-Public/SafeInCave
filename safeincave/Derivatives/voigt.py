# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tensorial-Voigt packing shared by the ISV-coupled consistent tangents."""

from __future__ import annotations
import torch as to

from .base import VOIGT_ROWS_I, VOIGT_ROWS_J


def pack_H_voigt(Q: to.Tensor, P: to.Tensor) -> to.Tensor:
    """
    Pack the rank-one tensor ``Q (outer) P`` into a 6×6 tensorial-Voigt matrix.

    ``H[k, l] = mult[l] * Q[rows[k]] * P[rows[l]]`` where `rows` selects the
    six independent components ``[xx, yy, zz, xy, xz, yz]`` and the three
    shear *columns* carry a factor of 2 — the multiplicity implicit in
    tensorial Voigt storage of a symmetric tensor.

    Parameters
    ----------
    Q : torch.Tensor
        Sensitivity of the strain rate to the internal variable, shape (N, 3, 3).
    P : torch.Tensor
        Sensitivity of the residue to stress, shape (N, 3, 3).

    Returns
    -------
    torch.Tensor
        `H`, shape (N, 6, 6).

    Notes
    -----
    Replaces the three byte-identical ``compute_H`` / ``_compute_H`` methods
    that used to live in `ViscoplasticDesai`, `MunsonDawsonCreep` and
    `ModifiedCamClayViscoplastic`.
    """
    q = Q[:, VOIGT_ROWS_I, VOIGT_ROWS_J]  # (N, 6)
    p = P[:, VOIGT_ROWS_I, VOIGT_ROWS_J]  # (N, 6)

    multiplicity = to.tensor([1.0, 1.0, 1.0, 2.0, 2.0, 2.0], dtype=P.dtype)
    return q[:, :, None] * (p * multiplicity)[:, None, :]
