# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
from .NonElasticElement import NonElasticElement


def _rate_kernel(xp: Any, stress: Any, A: Any) -> Any:
    """Namespace-generic linear-dashpot rate kernel."""
    return A[:, None, None] * stress


class LinearDashpot(NonElasticElement):
    """
    Linear viscous element (dashpot).

    Parameters
    ----------
    A : torch.Tensor
        Viscosity coefficient per element, shape (N,).
    name : str, optional
        Element name, by default "dashpot".
    derivative_method : str or DerivativeEvaluator, optional
        Derivative backend; see :class:`NonElasticElement`.
    """

    def __init__(
        self,
        A: to.Tensor,
        name: str = "dashpot",
        derivative_method: Any = None,
    ) -> None:
        super().__init__(A.shape[0], derivative_method=derivative_method)
        self.A = A
        self.name = name

    def rate_fn(self, xp: Any, stress: Any, phi1: float, Temp: Any) -> Any:
        """Namespace-generic strain-rate kernel (see :class:`NonElasticElement`)."""
        return _rate_kernel(xp, stress, self._cast(xp, self.A))

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        return_eps_ne: bool = False,
    ) -> None:
        """
        Compute linear viscous strain rate.

        Parameters
        ----------
        stress_vec : torch.Tensor
            Stress tensor per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (unused).
        Temp : torch.Tensor
            Temperature per element (unused).
        return_eps_ne : bool, default=False
            If True, return the rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.
        """
        eps_rate = _rate_kernel(to, stress_vec, self.A)
        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
