# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from .NonElasticElement import NonElasticElement


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

    Attributes
    ----------
    R : float
        Gas constant used (8.32).
    A, Q, d : torch.Tensor
        Material parameters, shape (N,).
    """

    def __init__(self, A: to.Tensor, d: to.Tensor, Q: to.Tensor, name: str = "creep") -> None:
        super().__init__(A.shape[0])
        self.R = 8.32
        self.Q = Q
        self.A = A
        self.d = d
        self.name = name

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
        s_xx = stress_vec[:, 0, 0]
        s_yy = stress_vec[:, 1, 1]
        s_zz = stress_vec[:, 2, 2]

        sigma_mean = (s_xx + s_yy + s_zz) / 3
        dev = stress_vec.clone()
        dev[:, 0, 0] = s_xx - sigma_mean
        dev[:, 1, 1] = s_yy - sigma_mean
        dev[:, 2, 2] = s_zz - sigma_mean

        A_bar = (self.A / self.d**3 / Temp) * to.exp(-self.Q / self.R / Temp)
        eps_rate = A_bar[:, None, None] * dev
        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
