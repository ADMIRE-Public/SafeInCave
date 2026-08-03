# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
from ...Derivatives import deviator, split_sym3
from .NonElasticElement import NonElasticElement


def _rate_kernel(
    xp: Any, stress: Any, Temp: Any, A: Any, d: Any, Q: Any, R: float
) -> Any:
    """
    Namespace-generic pressure-solution-creep rate kernel.

    Pure and out-of-place so it can be traced by both AD backends.
    """
    s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = split_sym3(stress)

    dev = deviator(xp, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

    A_bar = (A / d**3 / Temp) * xp.exp(-Q / R / Temp)
    return A_bar[:, None, None] * dev


class PressureSolutionCreep(NonElasticElement):
    """
    Pressure solution creep: :math:`\\dot\\varepsilon_{ne}
    = A/(Td^3)\\,\\exp(-Q/(RT))\\,\\mathbf{s}`.

    Parameters
    ----------
    A : torch.Tensor
        Pre-exponential factor per element, shape (N,).
    d : torch.Tensor
        Grain size (diameter), shape (N,).
    Q : torch.Tensor
        Activation energy per element, shape (N,).
    name : str, optional
        Element name, by default "creep".
    derivative_method : str or DerivativeEvaluator, optional
        Derivative backend; see :class:`NonElasticElement`.

    Attributes
    ----------
    R : float
        Gas constant used (8.32).
    A, Q, d : torch.Tensor
        Material parameters, shape (N,).
    """

    def __init__(
        self,
        A: to.Tensor,
        d: to.Tensor,
        Q: to.Tensor,
        name: str = "creep",
        derivative_method: Any = None,
    ) -> None:
        super().__init__(A.shape[0], derivative_method=derivative_method)
        self.R = 8.32
        self.Q = Q
        self.A = A
        self.d = d
        self.name = name

    def rate_fn(self, xp: Any, stress: Any, phi1: float, Temp: Any) -> Any:
        """Namespace-generic strain-rate kernel (see :class:`NonElasticElement`)."""
        A, d, Q = self._cast(xp, self.A, self.d, self.Q)
        return _rate_kernel(xp, stress, Temp, A, d, Q, self.R)

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        return_eps_ne: bool = False,
    ) -> None:
        """
        Compute creep strain rate from current stress.

        Parameters
        ----------
        stress_vec : torch.Tensor
            Stress tensor per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (unused here).
        Temp : torch.Tensor
            Temperature per element (N,) or broadcastable.
        return_eps_ne : bool, default=False
            If True, return the rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.
        """
        eps_rate = _rate_kernel(
            to, stress_vec, Temp, self.A, self.d, self.Q, self.R
        )
        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
