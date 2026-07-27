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

    def _jvp_batch(
        self, fn: Callable, x: to.Tensor, directions: to.Tensor
    ) -> tuple:
        """
        Batched forward-mode JVPs along multiple tangent directions.

        Parameters
        ----------
        fn : callable
            ``fn(xp, x)`` returning a tensor or a tuple of tensors.
        x : torch.Tensor
            Evaluation point, shape (N, 3, 3).
        directions : torch.Tensor
            Batch of tangent directions, shape (B, N, 3, 3) where B is the
            number of directions.

        Returns
        -------
        tuple
            ``(value, tangents)`` where `value` matches the output of `fn(xp, x)`
            and `tangents` has shape (B, ...) stacking the derivative along each
            direction.

        Notes
        -----
        Default implementation loops over directions and stacks results. Backends
        should override this for efficiency (e.g. with vmap).
        """
        results = [self._jvp(fn, x, directions[i]) for i in range(directions.shape[0])]
        values, derivs = zip(*results)
        return values[0], to.stack(derivs)

    def d_tensor_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-2,
    ) -> to.Tensor:
        n_elems = stress.shape[0]
        E = to.zeros((n_elems, 6, 6), dtype=to.float64)

        # Build batch of 6 tangent directions (one per STRESS_BASIS entry)
        directions = to.zeros((6, n_elems, 3, 3), dtype=stress.dtype, device=stress.device)
        for idx, (i, j, k, phi) in enumerate(STRESS_BASIS):
            directions[idx, :, i, j] = 1.0

        # Compute all 6 JVPs in one batched call
        _, d_fn_batch = self._jvp_batch(fn, stress, directions)

        # Scatter results into Jacobian columns
        for idx, (i, j, k, phi) in enumerate(STRESS_BASIS):
            E[:, :, k] = phi * d_fn_batch[idx, :, VOIGT_ROWS_I, VOIGT_ROWS_J]
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

        # Build batch of 6 tangent directions (one per STRESS_BASIS entry)
        directions = to.zeros((6, stress.shape[0], 3, 3), dtype=stress.dtype, device=stress.device)
        for idx, (i, j, _, _) in enumerate(STRESS_BASIS):
            directions[idx, :, i, j] = 1.0

        # Compute all 6 JVPs in one batched call
        _, d_fn_batch = self._jvp_batch(fn, stress, directions)

        # Scatter results into symmetric matrix
        for idx, (i, j, _, _) in enumerate(STRESS_BASIS):
            P[:, i, j] = d_fn_batch[idx]
            P[:, j, i] = d_fn_batch[idx]
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
