# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from .NonElasticElement import NonElasticElement


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

    Attributes
    ----------
    R : float
        Gas constant used (8.32).
    A, Q, n : torch.Tensor
        Material parameters, shape (N,).
    """

    def __init__(self, A: to.Tensor, Q: to.Tensor, n: to.Tensor, name: str = "creep"):
        super().__init__(A.shape[0])
        self.R = 8.32
        self.Q = Q
        self.A = A
        self.n = n
        self.name = name

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        return_eps_ne: bool = False,
    ):
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
        s_xx = stress_vec[:, 0, 0]
        s_yy = stress_vec[:, 1, 1]
        s_zz = stress_vec[:, 2, 2]
        s_xy = stress_vec[:, 0, 1]
        s_xz = stress_vec[:, 0, 2]
        s_yz = stress_vec[:, 1, 2]

        sigma_mean = (s_xx + s_yy + s_zz) / 3
        dev = stress_vec.clone()
        dev[:, 0, 0] = s_xx - sigma_mean
        dev[:, 1, 1] = s_yy - sigma_mean
        dev[:, 2, 2] = s_zz - sigma_mean

        q_vm = to.sqrt(
            0.5
            * (
                (s_xx - s_yy) ** 2
                + (s_xx - s_zz) ** 2
                + (s_yy - s_zz) ** 2
                + 6 * (s_xy**2 + s_xz**2 + s_yz**2)
            )
        )

        A_bar = self.A * to.exp(-self.Q / self.R / Temp) * q_vm ** (self.n - 1)
        eps_rate = A_bar[:, None, None] * dev
        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
