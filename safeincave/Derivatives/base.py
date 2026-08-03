# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Abstract interface shared by every derivative backend."""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any, Callable, Optional
import torch as to


#: Perturbation basis used by every stress derivative in the package:
#: ``(i, j, voigt_column, multiplicity)``.  Only the upper triangle is
#: perturbed — ``stress[:, j, i]`` is left untouched — and the three shear
#: directions carry a multiplicity of 2 to account for the symmetric
#: counterpart in tensorial Voigt storage.
STRESS_BASIS: tuple = (
    (0, 0, 0, 1.0),
    (1, 1, 1, 1.0),
    (2, 2, 2, 1.0),
    (0, 1, 3, 2.0),
    (0, 2, 4, 2.0),
    (1, 2, 5, 2.0),
)

#: Row indices selecting the six independent components of a symmetric
#: (N, 3, 3) output tensor, in tensorial Voigt order.
VOIGT_ROWS_I: list = [0, 1, 2, 0, 0, 1]
VOIGT_ROWS_J: list = [0, 1, 2, 1, 2, 2]


class DerivativeEvaluator(ABC):
    """
    Strategy for differentiating constitutive kernels.

    Every callable passed to these methods is *namespace-generic*: it takes
    the array module as its first argument, so the same kernel can be
    evaluated with `torch` (finite differences, torch AD) or `jax.numpy`
    (JAX AD).  See :mod:`safeincave.Derivatives.namespace`.

    Notes
    -----
    All three primitives differentiate along the **unsymmetrised** basis
    directions listed in :data:`STRESS_BASIS`.  This is deliberate: the
    finite-difference code these backends replace perturbs `stress[:, i, j]`
    alone and scales shear columns by 2, and the rate laws are not symmetric
    under an `(i, j) <-> (j, i)` swap.  Differentiating a symmetrised
    6-vector parameterisation instead would give a *different* matrix, not a
    more accurate one, and would break the existing consistent tangent.
    """

    name: str = "base"

    @abstractmethod
    def d_tensor_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-2,
    ) -> to.Tensor:
        """
        Jacobian of a tensor-valued kernel with respect to stress.

        Parameters
        ----------
        fn : callable
            ``fn(xp, stress) -> (N, 3, 3)``.
        stress : torch.Tensor
            Evaluation point, shape (N, 3, 3).
        step : float, default 1e-2
            Perturbation size, in Pa.  Ignored by the AD backends.

        Returns
        -------
        torch.Tensor
            Shape (N, 6, 6), tensorial Voigt ordering
            ``[xx, yy, zz, xy, xz, yz]``.
        """

    @abstractmethod
    def d_scalar_d_stress(
        self,
        fn: Callable,
        stress: to.Tensor,
        *,
        step: float = 1e-1,
        f0: Optional[to.Tensor] = None,
    ) -> to.Tensor:
        """
        Gradient of a scalar-valued kernel with respect to stress.

        Parameters
        ----------
        fn : callable
            ``fn(xp, stress) -> (N,)``.
        stress : torch.Tensor
            Evaluation point, shape (N, 3, 3).
        step : float, default 1e-1
            Perturbation size, in Pa.  Ignored by the AD backends.
        f0 : torch.Tensor, optional
            Value of `fn` at `stress`, if already known.  Lets the
            finite-difference backend skip the reference evaluation.

        Returns
        -------
        torch.Tensor
            Shape (N, 3, 3), symmetric (the upper triangle is mirrored).
        """

    @abstractmethod
    def d_d_isv(
        self,
        fn: Callable,
        x: to.Tensor,
        *,
        step: Any,
        f0: Any = None,
    ) -> Any:
        """
        Derivative with respect to a scalar internal state variable.

        Parameters
        ----------
        fn : callable
            ``fn(xp, x)`` returning a tensor or a tuple of tensors, each with
            leading dimension N.
        x : torch.Tensor
            Evaluation point, shape (N,).
        step : torch.Tensor or float
            Per-element perturbation, shape (N,).  Ignored by the AD backends.
        f0 : optional
            Value of `fn` at `x`, with the same structure as its output, if
            already known.

        Returns
        -------
        Same structure as the output of `fn`.
        """

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"{type(self).__name__}(name={self.name!r})"
