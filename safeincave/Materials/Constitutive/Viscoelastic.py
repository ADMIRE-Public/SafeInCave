# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from ..Utils import dotdot_torch
from .NonElasticElement import NonElasticElement


class Viscoelastic(NonElasticElement):
    """
    Kelvin–Voigt-type viscoelastic element.

    Parameters
    ----------
    eta : torch.Tensor
        Viscosity parameter per element, shape (N,).
    E : torch.Tensor
        Young's modulus per element, shape (N,).
    nu : torch.Tensor
        Poisson's ratio per element, shape (N,).
    name : str, optional
        Element name, by default "kelvin_voigt".

    Attributes
    ----------
    C1 : torch.Tensor
        Elastic stiffness in Voigt form, shape (N, 6, 6).
    eta, E, nu : torch.Tensor
        Material parameters, shape (N,).
    """

    def __init__(
        self, eta: to.Tensor, E: to.Tensor, nu: to.Tensor, name: str = "kelvin_voigt"
    ) -> None:
        super().__init__(E.shape[0])
        self.eta = eta
        self.E = E
        self.nu = nu
        self.name = name

        # Assemble C1 tensor (n_elems, 6, 6)
        self.C1 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        a0 = self.E / ((1 + self.nu) * (1 - 2 * self.nu))
        self.C1[:, 0, 0] = a0 * (1 - self.nu)
        self.C1[:, 1, 1] = a0 * (1 - self.nu)
        self.C1[:, 2, 2] = a0 * (1 - self.nu)
        self.C1[:, 3, 3] = a0 * (1 - 2 * self.nu)
        self.C1[:, 4, 4] = a0 * (1 - 2 * self.nu)
        self.C1[:, 5, 5] = a0 * (1 - 2 * self.nu)
        self.C1[:, 0, 1] = self.C1[:, 1, 0] = self.C1[:, 0, 2] = self.C1[
            :, 2, 0
        ] = self.C1[:, 2, 1] = self.C1[:, 1, 2] = a0 * self.nu

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        return_eps_ne: bool = False,
    ) -> None:
        """
        Compute viscoelastic strain rate (Kelvin–Voigt form).

        Parameters
        ----------
        stress_vec : torch.Tensor
            Stress tensor per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (dt*theta).
        Temp : torch.Tensor
            Temperature per element (unused here).
        return_eps_ne : bool, default=False
            If True, return the rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.
        """
        eps_ne_rate = dotdot_torch(
            self.G,
            stress_vec
            - dotdot_torch(self.C1, self.eps_ne_old + phi1 * self.eps_ne_rate_old),
        )
        if return_eps_ne:
            return eps_ne_rate.clone()
        else:
            self.eps_ne_rate = eps_ne_rate.clone()

    def compute_E(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> to.Tensor:
        """
        Closed-form 6×6 operator for viscoelasticity:
        `E = (eta*I + phi2*C1)^{-1}`.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3). (Unused here.)
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element (unused).

        Returns
        -------
        torch.Tensor
            `E` with shape (N, 6, 6).
        """
        phi2 = dt * (1 - theta)
        I_3x3 = to.eye(6, dtype=to.float64).unsqueeze(0).repeat(self.n_elems, 1, 1)
        E = to.linalg.inv(self.eta[:, None, None] * I_3x3 + phi2 * self.C1)
        return E
