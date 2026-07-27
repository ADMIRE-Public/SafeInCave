# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Finite-difference derivative backend — the package default."""

from __future__ import annotations
from typing import Any, Callable, Optional
import torch as to

from .base import DerivativeEvaluator, STRESS_BASIS, VOIGT_ROWS_I, VOIGT_ROWS_J
from .namespace import broadcast_step, tree_map2


class FiniteDifferenceEvaluator(DerivativeEvaluator):
    """
    Numerical differentiation by finite differences.

    Central differences for :meth:`d_tensor_d_stress` (the consistent-tangent
    operator `E`), forward differences elsewhere.  The default step sizes are
    the ones the constitutive models have always used, so this backend
    reproduces the historical numerics exactly.
    """

    name = "finite_difference"

    def d_tensor_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-2,
    ) -> to.Tensor:
        n_elems = stress.shape[0]
        E = to.zeros((n_elems, 6, 6), dtype=to.float64)
        stress_eps = stress.clone()

        for i, j, k, phi in STRESS_BASIS:
            stress_eps[:, i, j] += step
            f_plus = fn(to, stress_eps)
            stress_eps[:, i, j] -= 2 * step
            f_minus = fn(to, stress_eps)
            stress_eps[:, i, j] += step

            E[:, :, k] = (
                phi
                * (
                    f_plus[:, VOIGT_ROWS_I, VOIGT_ROWS_J]
                    - f_minus[:, VOIGT_ROWS_I, VOIGT_ROWS_J]
                )
                / (2 * step)
            )
        return E

    def d_scalar_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-1,
        f0: Optional[to.Tensor] = None,
    ) -> to.Tensor:
        if f0 is None:
            f0 = fn(to, stress)

        P = to.zeros_like(stress)
        stress_eps = stress.clone()
        for i, j, _, _ in STRESS_BASIS:
            stress_eps[:, i, j] += step
            f_plus = fn(to, stress_eps)
            P[:, i, j] = (f_plus - f0) / step
            P[:, j, i] = P[:, i, j]
            stress_eps[:, i, j] -= step
        return P

    def d_d_isv(
        self,
        fn: Callable,
        x: to.Tensor,
        *,
        step: Any,
        f0: Any = None,
    ) -> Any:
        if f0 is None:
            f0 = fn(to, x)
        f1 = fn(to, x + step)
        return tree_map2(lambda a, b: (a - b) / broadcast_step(step, a.dim()), f1, f0)
