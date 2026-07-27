# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Automatic differentiation backend built on JAX.

`jax` is an **optional** dependency (``pip install safeincave[jax]``); it is
imported lazily so that a jax-less environment is unaffected by importing
this module.  Selecting ``derivative_method="jax_ad"`` without jax installed
raises a clear ImportError.

Only mechanisms that provide a namespace-generic rate kernel — i.e. that
override :meth:`~safeincave.Materials.Constitutive.NonElasticElement.NonElasticElement.rate_fn`
— can use this backend; the others raise NotImplementedError telling you so.
"""

from __future__ import annotations
from typing import Any, Callable
import torch as to

from .autodiff import ADEvaluatorBase
from .namespace import to_namespace, from_namespace, tree_map

_JAX_MODULES: tuple | None = None


def _import_jax() -> tuple:
    """Import jax on first use and enable float64 (jax defaults to float32)."""
    global _JAX_MODULES
    if _JAX_MODULES is None:
        try:
            import jax
            import jax.numpy as jnp
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise ImportError(
                "derivative_method='jax_ad' requires JAX. Install it with "
                "`pip install safeincave[jax]` (jax[cpu]==0.4.35, pinned by "
                "the numpy<2 constraint of dolfinx 0.9.0)."
            ) from exc
        jax.config.update("jax_enable_x64", True)
        _JAX_MODULES = (jax, jnp)
    return _JAX_MODULES


def jax_is_available() -> bool:
    """True when the JAX backend can be used in this environment."""
    try:
        _import_jax()
    except ImportError:
        return False
    return True


class JaxADEvaluator(ADEvaluatorBase):
    """Exact derivatives via `jax.jvp`, with torch tensors in and out."""

    name = "jax_ad"

    def _jvp(self, fn: Callable, x: to.Tensor, v: to.Tensor) -> tuple:
        jax, jnp = _import_jax()

        x_jax = to_namespace(jnp, x)
        v_jax = to_namespace(jnp, v)

        value, derivative = jax.jvp(lambda a: fn(jnp, a), (x_jax,), (v_jax,))
        return tree_map(from_namespace, value), tree_map(from_namespace, derivative)

    def _jvp_batch(
        self, fn: Callable, x: to.Tensor, directions: to.Tensor
    ) -> tuple:
        """Batched JVP using vmap over the tangent directions, converting once."""
        jax, jnp = _import_jax()

        x_jax = to_namespace(jnp, x)
        directions_jax = to_namespace(jnp, directions)

        def single_jvp(v):
            return jax.jvp(lambda a: fn(jnp, a), (x_jax,), (v,))

        values, tangents = jax.vmap(single_jvp)(directions_jax)
        return tree_map(from_namespace, values[0]), tree_map(from_namespace, tangents)

    def to_jax(self, *arrays: Any) -> Any:
        """Convert torch tensors into jnp arrays (convenience for kernels)."""
        _, jnp = _import_jax()
        return to_namespace(jnp, *arrays)
