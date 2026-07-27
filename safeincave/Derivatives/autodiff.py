# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Shared machinery for the automatic-differentiation backends.

Both AD backends express every derivative as a small number of
**forward-mode directional derivatives** (JVPs) along the perturbation
directions in :data:`~safeincave.Derivatives.base.STRESS_BASIS`.  Using JVPs
rather than a full Jacobian keeps the AD path an exact `step -> 0` limit of
the finite-difference path it replaces, and costs six sweeps instead of the
twelve function evaluations central differences need.

Subclasses only implement :meth:`ADEvaluatorBase._jvp`.
"""

from __future__ import annotations
from abc import abstractmethod
from typing import Any, Callable, Optional
import torch as to

from .base import DerivativeEvaluator, STRESS_BASIS, VOIGT_ROWS_I, VOIGT_ROWS_J


class ADEvaluatorBase(DerivativeEvaluator):
    """Backend-agnostic assembly of the three primitives from JVPs."""

    name = "autodiff"

    @abstractmethod
    def _jvp(self, fn: Callable, x: to.Tensor, v: to.Tensor) -> tuple:
        """
        Forward-mode directional derivative of `fn` at `x` along `v`.

        Parameters
        ----------
        fn : callable
            ``fn(xp, x)`` returning a tensor or a tuple of tensors.
        x, v : torch.Tensor
            Evaluation point and tangent direction, same shape.

        Returns
        -------
        tuple
            ``(value, derivative)`` as torch tensors, matching the structure
            of the output of `fn`.
        """

    def d_tensor_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-2,
    ) -> to.Tensor:
        n_elems = stress.shape[0]
        E = to.zeros((n_elems, 6, 6), dtype=to.float64)

        for i, j, k, phi in STRESS_BASIS:
            direction = to.zeros_like(stress)
            direction[:, i, j] = 1.0
            _, d_fn = self._jvp(fn, stress, direction)
            E[:, :, k] = phi * d_fn[:, VOIGT_ROWS_I, VOIGT_ROWS_J]
        return E

    def d_scalar_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-1,
        f0: Optional[to.Tensor] = None,
    ) -> to.Tensor:
        P = to.zeros_like(stress)

        for i, j, _, _ in STRESS_BASIS:
            direction = to.zeros_like(stress)
            direction[:, i, j] = 1.0
            _, d_fn = self._jvp(fn, stress, direction)
            P[:, i, j] = d_fn
            P[:, j, i] = d_fn
        return P

    def d_d_isv(
        self,
        fn: Callable,
        x: to.Tensor,
        *,
        step: Any,
        f0: Any = None,
    ) -> Any:
        _, d_fn = self._jvp(fn, x, to.ones_like(x))
        return d_fn
