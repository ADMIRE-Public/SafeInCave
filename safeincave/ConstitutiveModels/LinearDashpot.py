# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from .NonElasticElement import NonElasticElement


class LinearDashpot(NonElasticElement):
    """
    Linear viscous element (dashpot).

    Parameters
    ----------
    A : torch.Tensor
        Viscosity coefficient per element, shape (N,).
    name : str, optional
        Element name, by default "dashpot".
    """

    def __init__(self, A: to.Tensor, name: str = "dashpot"):
        super().__init__(A.shape[0])
        self.A = A
        self.name = name

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        return_eps_ne: bool = False,
    ):
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
        eps_rate = self.A[:, None, None] * stress_vec
        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
