# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
from ...Derivatives import deviator, split_sym3
from .NonElasticElement import NonElasticElement


def _rate_kernel(
    xp: Any, stress: Any, Temp: Any, A: Any, Q: Any, n: Any, R: float
) -> Any:
    """
    Namespace-generic dislocation-creep rate kernel.

    Pure and out-of-place so it can be traced by both AD backends; see
    :mod:`safeincave.Derivatives.namespace`.
    """
    s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = split_sym3(stress)

    dev = deviator(xp, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

    q_vm = xp.sqrt(
        0.5
        * (
            (s_xx - s_yy) ** 2
            + (s_xx - s_zz) ** 2
            + (s_yy - s_zz) ** 2
            + 6 * (s_xy**2 + s_xz**2 + s_yz**2)
        )
    )

    A_bar = A * xp.exp(-Q / R / Temp) * q_vm ** (n - 1)
    return A_bar[:, None, None] * dev


class DislocationCreep(NonElasticElement):
    """
    Power-law dislocation creep: :math:`\\dot\\varepsilon_{ne}
    = A\\,\\exp(-Q/(RT))\\,\\q^{n-1}\\,\\mathbf{s}`.

    Parameters
    ----------
    A : torch.Tensor
        Pre-exponential factor per element, shape (N,).
    Q : torch.Tensor
        Activation energy per element, shape (N,).
    n : torch.Tensor
        Stress exponent per element, shape (N,).
    name : str, optional
        Element name, by default "creep".
    derivative_method : str or DerivativeEvaluator, optional
        Derivative backend; see :class:`NonElasticElement`.

    Attributes
    ----------
    R : float
        Gas constant used (8.32).
    A, Q, n : torch.Tensor
        Material parameters, shape (N,).
    """

    def __init__(
        self,
        A: to.Tensor,
        Q: to.Tensor,
        n: to.Tensor,
        name: str = "creep",
        derivative_method: Any = None,
    ) -> None:
        super().__init__(A.shape[0], derivative_method=derivative_method)
        self.R = 8.32
        self.Q = Q
        self.A = A
        self.n = n
        self.name = name

    def rate_fn(self, xp: Any, stress: Any, phi1: float, Temp: Any) -> Any:
        """Namespace-generic strain-rate kernel (see :class:`NonElasticElement`)."""
        A, Q, n = self._cast(xp, self.A, self.Q, self.n)
        return _rate_kernel(xp, stress, Temp, A, Q, n, self.R)

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
            to, stress_vec, Temp, self.A, self.Q, self.n, self.R
        )
        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
