# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
import numpy as np
from ..Utils import dotdot_torch, MPa
from .NonElasticElement import NonElasticElement


class ViscoplasticDesai(NonElasticElement):
    """
    Viscoplastic model of Desai-type with hardening (state) variable `alpha`.

    Parameters
    ----------
    mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0 : torch.Tensor
        Model parameters per element, shape (N,).
    name : str, optional
        Element name, by default "desai".

    Attributes
    ----------
    alpha : torch.Tensor
        Current hardening variable per element, shape (N,).
    Fvp : torch.Tensor
        Current value of the yield function per element, shape (N,).
    qsi, qsi_old : torch.Tensor
        Accumulated viscoplastic strain measure and its previous value, shape (N,).
    P : torch.Tensor
        Sensitivity of the residue to stress, shape (N, 3, 3).
    r, h : torch.Tensor
        Residue and its derivative w.r.t. `alpha`, shapes (N,) and (N,).
    """

    def __init__(
        self,
        mu_1: to.Tensor,
        N_1: to.Tensor,
        a_1: to.Tensor,
        eta: to.Tensor,
        n: to.Tensor,
        beta_1: to.Tensor,
        beta: to.Tensor,
        m: to.Tensor,
        gamma: to.Tensor,
        sigma_t: to.Tensor,
        alpha_0: to.Tensor,
        name: str = "desai",
    ):
        super().__init__(mu_1.shape[0])
        self.name = name
        self.mu_1 = mu_1
        self.N_1 = N_1
        self.a_1 = a_1
        self.eta = eta
        self.n = n
        self.beta_1 = beta_1
        self.beta = beta
        self.m = m
        self.gamma = gamma
        self.sigma_t = sigma_t
        self.alpha_0 = alpha_0
        self.F_0 = 1.0
        self.n_elems = self.alpha_0.shape[0]
        self.alpha = self.alpha_0.clone()
        self.Fvp = to.zeros(self.n_elems, dtype=to.float64)
        self.qsi = to.zeros(self.n_elems, dtype=to.float64)
        self.qsi_old = to.zeros(self.n_elems, dtype=to.float64)

    def compute_residue(
        self, eps_rate: to.Tensor, alpha: to.Tensor, dt: float
    ) -> to.Tensor:
        """
        Residue of the implicit hardening equation.

        Parameters
        ----------
        eps_rate : torch.Tensor
            Current inelastic strain rate, shape (N, 3, 3).
        alpha : torch.Tensor
            Hardening variable, shape (N,).
        dt : float
            Time step.

        Returns
        -------
        torch.Tensor
            Residue per element, shape (N,).

        Notes
        -----
        Updates `self.qsi` internally based on `eps_rate`.
        """
        self.qsi = self.qsi_old + to.sum(eps_rate**2, axis=(-2, -1)) ** 0.5 * dt
        return alpha - self.a_1 / (
            ((self.a_1 / self.alpha_0) ** (1 / self.eta) + self.qsi) ** self.eta
        )

    def update_internal_variables(self) -> None:
        """
        Commit accumulated measure `qsi`.

        Returns
        -------
        None
        """
        self.qsi_old = self.qsi.clone()

    def increment_internal_variables(
        self, stress: to.Tensor, stress_k: to.Tensor, dt: float
    ) -> None:
        """
        Increment hardening variable `alpha` using linearization.

        Parameters
        ----------
        stress : torch.Tensor
            End-of-step stress, shape (N, 3, 3).
        stress_k : torch.Tensor
            Intermediate stress, shape (N, 3, 3).
        dt : float
            Time step.

        Returns
        -------
        None

        Side Effects
        ------------
        Updates `self.alpha`.
        """
        delta_alpha = (
            -(self.r + to.einsum("bij,bij->b", self.P, stress - stress_k)) / self.h
        )
        self.alpha += delta_alpha

    def compute_stress_invariants(
        self,
        s_xx: to.Tensor,
        s_yy: to.Tensor,
        s_zz: to.Tensor,
        s_xy: to.Tensor,
        s_xz: to.Tensor,
        s_yz: to.Tensor,
    ) -> tuple[
        to.Tensor,
        to.Tensor,
        to.Tensor,
        to.Tensor,
        to.Tensor,
        to.Tensor,
        to.Tensor,
        to.Tensor,
    ]:
        """
        Compute invariants (I1, I2, I3, J2, J3) and auxiliary quantities.

        Parameters
        ----------
        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz : torch.Tensor
            Normal and shear components (MPa-scaled as provided), each shape (N,).

        Returns
        -------
        I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 : tuple of torch.Tensor
            Invariants and helper arrays; `ind_J2_leq_0` are indices where `J2 <= 0`.
        """
        I1 = s_xx + s_yy + s_zz
        I2 = s_xx * s_yy + s_yy * s_zz + s_xx * s_zz - s_xy**2 - s_yz**2 - s_xz**2
        I3 = (
            s_xx * s_yy * s_zz
            + 2 * s_xy * s_yz * s_xz
            - s_zz * s_xy**2
            - s_xx * s_yz**2
            - s_yy * s_xz**2
        )
        J2 = (1 / 3) * I1**2 - I2
        J3 = (2 / 27) * I1**3 - (1 / 3) * I1 * I2 + I3
        Sr = -(J3 * np.sqrt(27)) / (2 * J2**1.5)

        # Check where J2 <= 0.0
        ind_J2_leq_0 = to.where(J2 <= 0.0)[0]

        # Sr will be nan if, J2=0.0. So, replace it by 0.0
        Sr[ind_J2_leq_0] = 0.0

        I1_star = I1 + self.sigma_t
        return I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0

    def extract_stress_components(self, stress: to.Tensor) -> to.Tensor:
        """
        Extract and scale stress components from a 3×3 tensor.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).

        Returns
        -------
        tuple[torch.Tensor, ...]
            `(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)`, each shape (N,).
        """
        stress_vec = -stress
        s_xx = stress_vec[:, 0, 0] / MPa
        s_yy = stress_vec[:, 1, 1] / MPa
        s_zz = stress_vec[:, 2, 2] / MPa
        s_xy = stress_vec[:, 0, 1] / MPa
        s_xz = stress_vec[:, 0, 2] / MPa
        s_yz = stress_vec[:, 1, 2] / MPa
        return s_xx, s_yy, s_zz, s_xy, s_xz, s_yz

    def compute_Fvp(self, alpha, I1, J2, Sr):
        """
        Compute the Desai viscoplastic yield function value `Fvp`.

        Parameters
        ----------
        alpha : torch.Tensor
            Hardening variable per element, shape (N,).
        I1 : torch.Tensor
            First invariant, shape (N,).
        J2 : torch.Tensor
            Second deviatoric invariant, shape (N,).
        Sr : torch.Tensor
            Lode-related parameter, shape (N,).

        Returns
        -------
        torch.Tensor
            `Fvp` per element, shape (N,).
        """
        F1 = alpha * I1**self.n - self.gamma * I1**2
        F2 = to.exp(self.beta_1 * I1) - self.beta * Sr
        Fvp = J2 + F1 * F2**self.m
        return Fvp

    def compute_initial_hardening(self, stress: to.Tensor, Fvp_0=0.0) -> None:
        """
        Initialize `alpha` from a target `Fvp_0` and the current stress state.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        Fvp_0 : float, default=0.0
            Target initial value for `Fvp`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.alpha_0`, `self.alpha`, and `self.Fvp`.
        """
        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
        I1, I2, I3, J2, J3, Sr, I1_star, _ = self.compute_stress_invariants(
            s_xx, s_yy, s_zz, s_xy, s_xz, s_yz
        )
        self.alpha_0 = self.gamma * I1_star ** (2 - self.n) + (
            Fvp_0 - J2
        ) * I1_star ** (-self.n) * (to.exp(self.beta_1 * I1_star) - self.beta * Sr) ** (
            -self.m
        )
        self.alpha = self.alpha_0.clone()

        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
        I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 = self.compute_stress_invariants(
            s_xx, s_yy, s_zz, s_xy, s_xz, s_yz
        )
        self.Fvp = self.compute_Fvp(self.alpha, I1_star, J2, Sr)

    def compute_eps_ne_rate(
        self,
        stress: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        alpha=None,
        return_eps_ne=False,
    ):
        """
        Compute viscoplastic strain rate and optionally return it.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (dt*theta).
        Temp : torch.Tensor
            Temperature per element (unused).
        alpha : torch.Tensor or None, optional
            Hardening variable override; if `None`, use `self.alpha`.
        return_eps_ne : bool, default=False
            If True, return rate; else store it.

        Returns
        -------
        None or torch.Tensor
            (N, 3, 3) if `return_eps_ne=True`, else `None`.

        Notes
        -----
        Also updates `self.Fvp` when `return_eps_ne=False`.
        """
        if alpha is None:
            alpha = self.alpha

        s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = self.extract_stress_components(stress)
        I1, I2, I3, J2, J3, Sr, I1_star, ind_J2_leq_0 = self.compute_stress_invariants(
            s_xx, s_yy, s_zz, s_xy, s_xz, s_yz
        )

        # Compute yield function
        Fvp = self.compute_Fvp(alpha, I1_star, J2, Sr)
        if not return_eps_ne:
            self.Fvp = Fvp.clone()

        # Compute flow direction, i.e. d(Fvp)/d(stress)
        F1 = -alpha * I1_star**self.n + self.gamma * I1_star**2
        F2 = to.exp(self.beta_1 * I1_star) - self.beta * Sr
        dF1_dI1 = 2 * self.gamma * I1_star - self.n * alpha * I1_star ** (self.n - 1)
        dF2m_dI1 = (
            self.beta_1 * self.m * to.exp(self.beta_1 * I1_star) * F2 ** (self.m - 1)
        )
        dF_dI1 = -(dF1_dI1 * F2**self.m + F1 * dF2m_dI1)

        dF2_dJ2 = -(3 * self.beta * J3 * 27**0.5) / (4 * J2 ** (5 / 2))
        dF_dJ2 = 1 - F1 * self.m * F2 ** (self.m - 1) * dF2_dJ2
        dF_dJ3 = (
            -self.m * F1 * self.beta * np.sqrt(27) * F2 ** (self.m - 1) / (2 * J2**1.5)
        )

        dI1_dSxx = 1.0
        dI1_dSyy = 1.0
        dI1_dSzz = 1.0
        dI1_dSxy = 0.0
        dI1_dSxz = 0.0
        dI1_dSyz = 0.0

        dI2_dSxx = s_yy + s_zz
        dI2_dSyy = s_xx + s_zz
        dI2_dSzz = s_xx + s_yy
        dI2_dSxy = -2 * s_xy
        dI2_dSxz = -2 * s_xz
        dI2_dSyz = -2 * s_yz

        dI3_dSxx = s_yy * s_zz - s_yz**2
        dI3_dSyy = s_xx * s_zz - s_xz**2
        dI3_dSzz = s_xx * s_yy - s_xy**2
        dI3_dSxy = 2 * (s_xz * s_yz - s_zz * s_xy)
        dI3_dSxz = 2 * (s_xy * s_yz - s_yy * s_xz)
        dI3_dSyz = 2 * (s_xz * s_xy - s_xx * s_yz)

        dJ2_dI1 = (2 / 3) * I1
        dJ2_dI2 = -1.0

        dJ2_dSxx = dJ2_dI1 * dI1_dSxx + dJ2_dI2 * dI2_dSxx
        dJ2_dSyy = dJ2_dI1 * dI1_dSyy + dJ2_dI2 * dI2_dSyy
        dJ2_dSzz = dJ2_dI1 * dI1_dSzz + dJ2_dI2 * dI2_dSzz
        dJ2_dSxy = dJ2_dI1 * dI1_dSxy + dJ2_dI2 * dI2_dSxy
        dJ2_dSxz = dJ2_dI1 * dI1_dSxz + dJ2_dI2 * dI2_dSxz
        dJ2_dSyz = dJ2_dI1 * dI1_dSyz + dJ2_dI2 * dI2_dSyz

        dJ3_dI1 = (2 / 9) * I1**2 - (1 / 3) * I2
        dJ3_dI2 = -(1 / 3) * I1
        dJ3_dI3 = 1.0

        dJ3_dSxx = dJ3_dI1 * dI1_dSxx + dJ3_dI2 * dI2_dSxx + dJ3_dI3 * dI3_dSxx
        dJ3_dSyy = dJ3_dI1 * dI1_dSyy + dJ3_dI2 * dI2_dSyy + dJ3_dI3 * dI3_dSyy
        dJ3_dSzz = dJ3_dI1 * dI1_dSzz + dJ3_dI2 * dI2_dSzz + dJ3_dI3 * dI3_dSzz
        dJ3_dSxy = dJ3_dI1 * dI1_dSxy + dJ3_dI2 * dI2_dSxy + dJ3_dI3 * dI3_dSxy
        dJ3_dSxz = dJ3_dI1 * dI1_dSxz + dJ3_dI2 * dI2_dSxz + dJ3_dI3 * dI3_dSxz
        dJ3_dSyz = dJ3_dI1 * dI1_dSyz + dJ3_dI2 * dI2_dSyz + dJ3_dI3 * dI3_dSyz

        dQdS_00 = dF_dI1 * dI1_dSxx + dF_dJ2 * dJ2_dSxx + dF_dJ3 * dJ3_dSxx
        dQdS_11 = dF_dI1 * dI1_dSyy + dF_dJ2 * dJ2_dSyy + dF_dJ3 * dJ3_dSyy
        dQdS_22 = dF_dI1 * dI1_dSzz + dF_dJ2 * dJ2_dSzz + dF_dJ3 * dJ3_dSzz
        dQdS_01 = dF_dI1 * dI1_dSxy + dF_dJ2 * dJ2_dSxy + dF_dJ3 * dJ3_dSxy
        dQdS_02 = dF_dI1 * dI1_dSxz + dF_dJ2 * dJ2_dSxz + dF_dJ3 * dJ3_dSxz
        dQdS_12 = dF_dI1 * dI1_dSyz + dF_dJ2 * dJ2_dSyz + dF_dJ3 * dJ3_dSyz

        # Initialize viscoplastic direction
        dQdS = to.zeros_like(stress, dtype=to.float64)
        dQdS[:, 0, 0] = dQdS_00
        dQdS[:, 1, 1] = dQdS_11
        dQdS[:, 2, 2] = dQdS_22
        dQdS[:, 1, 0] = dQdS[:, 0, 1] = dQdS_01
        dQdS[:, 2, 0] = dQdS[:, 0, 2] = dQdS_02
        dQdS[:, 2, 1] = dQdS[:, 1, 2] = dQdS_12

        # Wherever J2=0, make viscoplasticity to be zero
        dQdS[ind_J2_leq_0, :, :] = 0.0

        # Calculate strain rate
        ramp_idx = to.where(Fvp > 0)[0]
        lmbda = to.zeros(self.n_elems, dtype=to.float64)
        if len(ramp_idx) != 0:
            lmbda[ramp_idx] = (
                self.mu_1[ramp_idx] * (Fvp[ramp_idx] / self.F_0) ** self.N_1[ramp_idx]
            )
        eps_vp_rate = -dQdS * lmbda[:, None, None]

        if return_eps_ne:
            return eps_vp_rate
        else:
            self.eps_ne_rate = eps_vp_rate

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> tuple[to.Tensor, to.Tensor]:
        """
        Compute `B` and `H/h` via perturbations of `alpha` and stress.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element (unused).

        Returns
        -------
        B : torch.Tensor
            Driving term, shape (N, 3, 3).
        H_over_h : torch.Tensor
            Linearization ratio, shape (N, 6, 6).

        Notes
        -----
        Uses finite differences to approximate sensitivities.
        """
        # EPSILON_ALPHA = 1e-7
        EPSILON_ALPHA = 0.0001 * self.alpha
        EPSILON_STRESS = 1e-1

        alpha_eps = self.alpha + EPSILON_ALPHA
        eps_ne_rate_eps = self.compute_eps_ne_rate(
            stress, dt * theta, Temp, alpha=alpha_eps, return_eps_ne=True
        )

        self.r = self.compute_residue(self.eps_ne_rate, self.alpha, dt)
        r_eps = self.compute_residue(eps_ne_rate_eps, alpha_eps, dt)
        self.h = (r_eps - self.r) / EPSILON_ALPHA
        Q = (eps_ne_rate_eps - self.eps_ne_rate) / EPSILON_ALPHA[:, None, None]
        B = (self.r / self.h)[:, None, None] * Q

        self.P = to.zeros_like(stress)
        stress_eps = stress.clone()
        for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
            stress_eps[:, i, j] += EPSILON_STRESS
            eps_ne_rate_eps = self.compute_eps_ne_rate(
                stress_eps, dt * theta, Temp, return_eps_ne=True
            )
            r_eps = self.compute_residue(eps_ne_rate_eps, self.alpha, dt)
            self.P[:, i, j] = (r_eps - self.r) / EPSILON_STRESS
            self.P[:, j, i] = self.P[:, i, j]
            stress_eps[:, i, j] -= EPSILON_STRESS

        H = self.compute_H(Q, self.P)
        H_over_h = H / self.h[:, None, None]

        return B, H_over_h

    def compute_H(self, Q: to.Tensor, P: to.Tensor) -> to.Tensor:
        """
        Build the 6×6 matrix `H` from tensors `Q` and `P`.

        Parameters
        ----------
        Q : torch.Tensor
            Sensitivity of rate to `alpha`, shape (N, 3, 3).
        P : torch.Tensor
            Sensitivity of residue to stress, shape (N, 3, 3).

        Returns
        -------
        torch.Tensor
            `H` with shape (N, 6, 6) in tensorial Voigt ordering.
        """
        n_elems, _, _ = P.shape
        H = to.zeros((n_elems, 6, 6), dtype=to.float64)
        H[:, 0, 0] = Q[:, 0, 0] * P[:, 0, 0]
        H[:, 0, 1] = Q[:, 0, 0] * P[:, 1, 1]
        H[:, 0, 2] = Q[:, 0, 0] * P[:, 2, 2]
        H[:, 0, 3] = 2 * Q[:, 0, 0] * P[:, 0, 1]
        H[:, 0, 4] = 2 * Q[:, 0, 0] * P[:, 0, 2]
        H[:, 0, 5] = 2 * Q[:, 0, 0] * P[:, 1, 2]

        H[:, 1, 0] = Q[:, 1, 1] * P[:, 0, 0]
        H[:, 1, 1] = Q[:, 1, 1] * P[:, 1, 1]
        H[:, 1, 2] = Q[:, 1, 1] * P[:, 2, 2]
        H[:, 1, 3] = 2 * Q[:, 1, 1] * P[:, 0, 1]
        H[:, 1, 4] = 2 * Q[:, 1, 1] * P[:, 0, 2]
        H[:, 1, 5] = 2 * Q[:, 1, 1] * P[:, 1, 2]

        H[:, 2, 0] = Q[:, 2, 2] * P[:, 0, 0]
        H[:, 2, 1] = Q[:, 2, 2] * P[:, 1, 1]
        H[:, 2, 2] = Q[:, 2, 2] * P[:, 2, 2]
        H[:, 2, 3] = 2 * Q[:, 2, 2] * P[:, 0, 1]
        H[:, 2, 4] = 2 * Q[:, 2, 2] * P[:, 0, 2]
        H[:, 2, 5] = 2 * Q[:, 2, 2] * P[:, 1, 2]

        H[:, 3, 0] = Q[:, 0, 1] * P[:, 0, 0]
        H[:, 3, 1] = Q[:, 0, 1] * P[:, 1, 1]
        H[:, 3, 2] = Q[:, 0, 1] * P[:, 2, 2]
        H[:, 3, 3] = 2 * Q[:, 0, 1] * P[:, 0, 1]
        H[:, 3, 4] = 2 * Q[:, 0, 1] * P[:, 0, 2]
        H[:, 3, 5] = 2 * Q[:, 0, 1] * P[:, 1, 2]

        H[:, 4, 0] = Q[:, 0, 2] * P[:, 0, 0]
        H[:, 4, 1] = Q[:, 0, 2] * P[:, 1, 1]
        H[:, 4, 2] = Q[:, 0, 2] * P[:, 2, 2]
        H[:, 4, 3] = 2 * Q[:, 0, 2] * P[:, 0, 1]
        H[:, 4, 4] = 2 * Q[:, 0, 2] * P[:, 0, 2]
        H[:, 4, 5] = 2 * Q[:, 0, 2] * P[:, 1, 2]

        H[:, 5, 0] = Q[:, 1, 2] * P[:, 0, 0]
        H[:, 5, 1] = Q[:, 1, 2] * P[:, 1, 1]
        H[:, 5, 2] = Q[:, 1, 2] * P[:, 2, 2]
        H[:, 5, 3] = 2 * Q[:, 1, 2] * P[:, 0, 1]
        H[:, 5, 4] = 2 * Q[:, 1, 2] * P[:, 0, 2]
        H[:, 5, 5] = 2 * Q[:, 1, 2] * P[:, 1, 2]
        return H
