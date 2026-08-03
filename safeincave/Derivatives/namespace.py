# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Array-namespace helpers shared by the derivative backends.

Constitutive rate laws are written once as *namespace-generic* kernels that
take the array module (`xp`) as their first argument:

    def _rate_kernel(xp, stress, ...):
        ...
        return build_sym3(xp, exx, eyy, ezz, exy, exz, eyz)

`xp` is `torch` for the finite-difference and torch-AD backends and
`jax.numpy` for the JAX-AD backend.  Only the small intersection of the two
APIs may be used inside a kernel: ``sqrt``, ``exp``, ``log10``, ``abs``,
``where``, ``clip``, ``stack``, ``sum``, ``ones_like``, ``zeros_like`` and
arithmetic.  In particular **no in-place assignment** (jnp arrays are
immutable) and **no data-dependent control flow** (breaks tracing and
silently zeroes gradients).
"""

from __future__ import annotations
from typing import Any
import numpy as np
import torch as to


def is_torch(xp: Any) -> bool:
    """True when `xp` is the PyTorch namespace."""
    return xp is to


def to_namespace(xp: Any, *arrays: Any) -> Any:
    """
    Convert torch tensors (or scalars) into the namespace `xp`.

    Returns a single value when called with one argument, otherwise a tuple.
    Scalars and non-tensor values pass through untouched.
    """
    converted = tuple(_convert_one(xp, a) for a in arrays)
    return converted[0] if len(converted) == 1 else converted


def _convert_one(xp: Any, a: Any) -> Any:
    if is_torch(xp) or not isinstance(a, to.Tensor):
        return a
    return xp.asarray(a.detach().cpu().numpy())


def from_namespace(a: Any) -> Any:
    """Convert an array produced by any namespace back into a torch tensor."""
    if isinstance(a, to.Tensor):
        return a
    return to.from_numpy(np.asarray(a, dtype=np.float64))


def build_sym3(
    xp: Any,
    xx: Any,
    yy: Any,
    zz: Any,
    xy: Any,
    xz: Any,
    yz: Any,
) -> Any:
    """
    Assemble a symmetric (N, 3, 3) tensor from its six independent components.

    Out-of-place (uses `stack`), so the result is differentiable under both
    torch forward-mode AD and JAX.
    """
    row0 = xp.stack([xx, xy, xz], -1)
    row1 = xp.stack([xy, yy, yz], -1)
    row2 = xp.stack([xz, yz, zz], -1)
    return xp.stack([row0, row1, row2], -2)


def split_sym3(stress: Any) -> tuple:
    """
    Split an (N, 3, 3) tensor into `(xx, yy, zz, xy, xz, yz)`.

    Reads only the upper triangle, matching how every rate law in the package
    consumes stress.
    """
    return (
        stress[:, 0, 0],
        stress[:, 1, 1],
        stress[:, 2, 2],
        stress[:, 0, 1],
        stress[:, 0, 2],
        stress[:, 1, 2],
    )


def deviator(xp: Any, xx: Any, yy: Any, zz: Any, xy: Any, xz: Any, yz: Any) -> Any:
    """Deviatoric part of a symmetric tensor given its six components."""
    mean = (xx + yy + zz) / 3.0
    return build_sym3(xp, xx - mean, yy - mean, zz - mean, xy, xz, yz)


def tree_map(func, obj: Any) -> Any:
    """Apply `func` to every array in a tensor / tuple / list structure."""
    if isinstance(obj, tuple):
        return tuple(tree_map(func, o) for o in obj)
    if isinstance(obj, list):
        return [tree_map(func, o) for o in obj]
    return func(obj)


def tree_map2(func, a: Any, b: Any) -> Any:
    """Apply `func` pairwise over two matching tensor / tuple / list structures."""
    if isinstance(a, tuple):
        return tuple(tree_map2(func, x, y) for x, y in zip(a, b))
    if isinstance(a, list):
        return [tree_map2(func, x, y) for x, y in zip(a, b)]
    return func(a, b)


def broadcast_step(step: Any, ndim: int) -> Any:
    """Right-pad a per-element step `(N,)` so it broadcasts against `ndim` axes."""
    if not isinstance(step, to.Tensor):
        return step
    while step.dim() < ndim:
        step = step[..., None]
    return step
