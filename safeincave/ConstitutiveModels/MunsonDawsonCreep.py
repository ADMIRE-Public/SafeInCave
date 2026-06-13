# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from .NonElasticElement import NonElasticElement


class MunsonDawsonCreep(NonElasticElement):
    """
    Munson–Dawson creep law (steady-state + transient) with internal variable zeta.

    Constitutive equations (formula sheet):
        epsdot_MD_ij = F * epsdot_ss * (3/2) s_ij / sigma
        zeta_dot     = (F - 1) * epsdot_ss
        epsdot_ss    = A * exp(-Q/(R T)) * sigma^n
        sigma        = sqrt(3 J2)
        eps_t*       = K0 * exp(c T) * (sigma/mu)^m
        Delta        = alpha_w + beta_w * log10(sigma/mu)
        F = exp( +Delta * (1 - zeta/eps_t*)^2 )   when zeta <= eps_t*   (hardening)
        F = exp( -delta * (1 - zeta/eps_t*)^2 )   when zeta  > eps_t*   (recovery)

    Numerical scheme
    ----------------
    zeta is treated as a true internal state variable, with the same
    consistent-tangent pattern used by `ViscoplasticDesai`.  Its update is
    linearised into the global Newton iteration via the residue equation

        r(zeta, sigma) = zeta - zeta_old - (F(zeta, sigma) - 1) * epsdot_ss(sigma) * dt

    Parameters
    ----------
    A, Q, n, K0, c, m, alpha_w, beta_w, delta, mu : torch.Tensor
        Per-element Munson-Dawson parameters, shape (N,).
    name : str, optional
        Element name, default "creep_munson_dawson".

    Attributes
    ----------
    zeta, zeta_old : torch.Tensor, shape (N,)
        Transient internal variable (current iterate and last-committed).
    F, _eps_t_star : torch.Tensor, shape (N,)
        Diagnostics cached by `compute_eps_ne_rate`.
    r, h, P : torch.Tensor
        Residue, dr/dzeta, and dr/dsigma populated by
        `compute_B_and_H_over_h` and consumed by
        `increment_internal_variables`.
    """

    # Optimal forward-difference step factor for double precision.
    _SQRT_FLOAT64_EPS = 1.4901161193847656e-8  # math.sqrt(2.220446e-16)

    def __init__(
        self,
        A: to.Tensor,
        Q: to.Tensor,
        n: to.Tensor,
        K0: to.Tensor,
        c: to.Tensor,
        m: to.Tensor,
        alpha_w: to.Tensor,
        beta_w: to.Tensor,
        delta: to.Tensor,
        mu: to.Tensor,
        name: str = "creep_munson_dawson",
    ):
        super().__init__(A.shape[0])
        self.name = name

        self.R = 8.32

        # Steady-state params
        self.A = A.to(dtype=to.float64)
        self.Q = Q.to(dtype=to.float64)
        self.n = n.to(dtype=to.float64)

        # Transient params
        self.K0 = K0.to(dtype=to.float64)
        self.c = c.to(dtype=to.float64)
        self.m = m.to(dtype=to.float64)
        self.alpha_w = alpha_w.to(dtype=to.float64)
        self.beta_w = beta_w.to(dtype=to.float64)
        self.delta = delta.to(dtype=to.float64)

        # Shear modulus (Pa)
        self.mu = mu.to(dtype=to.float64)

        # Internal variable zeta (starts at 0 – untransient-strained rock)
        self.zeta = to.zeros(self.n_elems, dtype=to.float64)
        self.zeta_old = self.zeta.clone()

        # Diagnostics populated by compute_eps_ne_rate
        self.F = to.ones(self.n_elems, dtype=to.float64)
        self._eps_t_star = to.ones(self.n_elems, dtype=to.float64)

        # Newton-coupling storage: filled by compute_B_and_H_over_h,
        # consumed by increment_internal_variables.
        self.r = to.zeros(self.n_elems, dtype=to.float64)
        self.h = to.ones(self.n_elems, dtype=to.float64)
        self.P = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.ind_h_small = to.tensor([], dtype=to.long)

    # ------------------------------------------------------------------ #
    # ISV lifecycle
    # ------------------------------------------------------------------ #

    def update_internal_variables(self) -> None:
        """Commit zeta at end of a converged time step."""
        self.zeta_old = self.zeta.clone()

    def increment_internal_variables(
        self, stress: to.Tensor, stress_k: to.Tensor, dt: float
    ) -> None:
        """
        Apply the linearised Newton correction for zeta using the stored
        residue and sensitivities from the most recent
        `compute_B_and_H_over_h`:

            delta_zeta = -(r + P : (stress - stress_k)) / h
        """
        delta_sigma = stress - stress_k
        delta_zeta = -(self.r + to.einsum("bij,bij->b", self.P, delta_sigma)) / self.h
        # Ill-conditioned elements: skip the update (P was zeroed, so the
        # contribution from delta_sigma vanishes, but guard the division too).
        if len(self.ind_h_small) > 0:
            delta_zeta[self.ind_h_small] = 0.0
        self.zeta = self.zeta + delta_zeta
        # zeta is physically non-negative.  Recovery (F < 1) can push it
        # back toward eps_t*, never below zero.
        self.zeta = to.clamp(self.zeta, min=0.0)

    # ------------------------------------------------------------------ #
    # Physics evaluator (shared by rate, residue, and FD probes)
    # ------------------------------------------------------------------ #

    def _compute_md_fields(
        self, stress_vec: to.Tensor, Temp: to.Tensor, zeta: to.Tensor
    ):
        """
        Compute all Munson-Dawson intermediate quantities for a given
        (stress, zeta) state.

        Returns
        -------
        s_dev : (N, 3, 3) deviatoric stress
        sigma_safe : (N,) von Mises stress (clamped at 1 Pa)
        epsdot_ss : (N,) steady-state scalar creep rate
        eps_t_star : (N,) transient threshold strain
        F : (N,) transient function value
        """
        s_xx = stress_vec[:, 0, 0]
        s_yy = stress_vec[:, 1, 1]
        s_zz = stress_vec[:, 2, 2]
        s_xy = stress_vec[:, 0, 1]
        s_xz = stress_vec[:, 0, 2]
        s_yz = stress_vec[:, 1, 2]

        sigma_mean = (s_xx + s_yy + s_zz) / 3.0
        s_dev = stress_vec.clone()
        s_dev[:, 0, 0] = s_xx - sigma_mean
        s_dev[:, 1, 1] = s_yy - sigma_mean
        s_dev[:, 2, 2] = s_zz - sigma_mean

        # sigma = sqrt(3 J2) = von Mises equivalent stress (Pa)
        sigma = to.sqrt(
            0.5
            * (
                (s_xx - s_yy) ** 2
                + (s_xx - s_zz) ** 2
                + (s_yy - s_zz) ** 2
                + 6.0 * (s_xy**2 + s_xz**2 + s_yz**2)
            )
        )

        # 1 Pa floor — purely numerical guard against 1/sigma and sigma^n at
        # zero-stress elements.  Does not affect any realistic cavern state.
        sigma_safe = to.clamp(sigma, min=1.0)
        mu_safe = to.clamp(self.mu, min=1.0)

        epsdot_ss = self.A * to.exp(-self.Q / (self.R * Temp)) * (sigma_safe**self.n)

        ratio = to.clamp(sigma_safe / mu_safe, min=1e-30)
        eps_t_star = self.K0 * to.exp(self.c * Temp) * (ratio**self.m)
        eps_t_star = to.clamp(eps_t_star, min=1e-50)

        Delta_cap = self.alpha_w + self.beta_w * to.log10(ratio)

        # Piecewise F with hardening / recovery branches.  The ±50 exponent
        # clamp guards against overflow during FD probes at pathological
        # stress states (exp(50) ≈ 5e21 — far beyond any realistic F).
        r_arg = 1.0 - (zeta / eps_t_star)
        r_arg2 = r_arg * r_arg

        F = to.ones_like(epsdot_ss)
        mask = zeta <= eps_t_star
        exp_arg_hard = to.clamp(Delta_cap[mask] * r_arg2[mask], min=-50.0, max=50.0)
        exp_arg_recov = to.clamp(
            -self.delta[~mask] * r_arg2[~mask], min=-50.0, max=50.0
        )
        F[mask] = to.exp(exp_arg_hard)
        F[~mask] = to.exp(exp_arg_recov)

        return s_dev, sigma_safe, epsdot_ss, eps_t_star, F

    def compute_residue(
        self, stress: to.Tensor, zeta: to.Tensor, Temp: to.Tensor, dt: float
    ) -> to.Tensor:
        """
        Backward-Euler residue for the zeta ODE:

            r(zeta, sigma) = zeta - zeta_old - (F(zeta, sigma) - 1) * epsdot_ss(sigma) * dt

        Returns
        -------
        torch.Tensor, shape (N,)
        """
        _, _, epsdot_ss, _, F = self._compute_md_fields(stress, Temp, zeta)
        return zeta - self.zeta_old - (F - 1.0) * epsdot_ss * dt

    # ------------------------------------------------------------------ #
    # Strain-rate evaluator
    # ------------------------------------------------------------------ #

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        zeta: to.Tensor = None,
        return_eps_ne: bool = False,
    ):
        """
        Munson-Dawson creep strain rate at (sigma, zeta).

        Parameters
        ----------
        stress_vec : torch.Tensor, shape (N, 3, 3)
        phi1 : float
            dt * theta — only present to match the base-class signature; the
            raw MD rate does not depend on it.
        Temp : torch.Tensor, shape (N,)
        zeta : torch.Tensor or None, optional
            Override zeta (used by the FD probe for Q = d eps_MD / d zeta).
            If None, uses `self.zeta`.
        return_eps_ne : bool, default False
            If True, return the rate without touching `self.eps_ne_rate`;
            required for FD probes.
        """
        if zeta is None:
            zeta = self.zeta

        s_dev, sigma_safe, epsdot_ss, eps_t_star, F = self._compute_md_fields(
            stress_vec, Temp, zeta
        )

        scalar_rate = F * epsdot_ss
        flow_dir = (1.5 / sigma_safe)[:, None, None] * s_dev
        eps_rate = flow_dir * scalar_rate[:, None, None]

        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
            # Cache diagnostics only on the "real" call (not FD probes,
            # which go through return_eps_ne=True).
            self._eps_t_star = eps_t_star
            self.F = F.clone()

    # ------------------------------------------------------------------ #
    # ISV-coupled consistent tangent
    # ------------------------------------------------------------------ #

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ):
        """
        Build the ISV-coupled consistent-tangent terms by finite differences
        (following ViscoplasticDesai.compute_B_and_H_over_h).

            r  = residue at the reference (sigma, zeta)
            h  = (r(sigma, zeta+eps) - r) / eps                (N,)
            Q  = (eps_MD(sigma, zeta+eps) - eps_MD) / eps      (N, 3, 3)
            P[i,j] = (r(sigma+eps*e_ij, zeta) - r) / eps       (N, 3, 3, symmetric)
            H   = Q (outer) P  packed in tensorial-Voigt 6×6
            H/h = (1/h) * H
            B   = (r / h) * Q

        r, h, P are stored on the instance for the subsequent
        `increment_internal_variables` call.  B and H/h are returned.
        """
        # Forward-difference scale for ζ: tied to the natural range eps_t*
        # so it stays meaningful whether ζ is 0 (start of run) or near ε_t*.
        _, _, _, eps_t_star_now, _ = self._compute_md_fields(stress, Temp, self.zeta)
        zeta_scale = to.clamp(to.abs(self.zeta) + eps_t_star_now, min=1e-30)
        eps_zeta = self._SQRT_FLOAT64_EPS * zeta_scale  # shape (N,)

        # Forward-difference scale for stress (Pa).  On ~1e7 Pa cavern
        # stresses this is a relative step of ~1e-8, which is the optimal
        # float64 FD scale.
        EPSILON_STRESS = 1e-1

        # --- r, h, Q via zeta-perturbation ------------------------------- #
        self.r = self.compute_residue(stress, self.zeta, Temp, dt)

        zeta_eps = self.zeta + eps_zeta
        r_zeta = self.compute_residue(stress, zeta_eps, Temp, dt)
        self.h = (r_zeta - self.r) / eps_zeta

        eps_rate_ref = self.compute_eps_ne_rate(
            stress, dt * theta, Temp, zeta=self.zeta, return_eps_ne=True
        )
        eps_rate_zeta = self.compute_eps_ne_rate(
            stress, dt * theta, Temp, zeta=zeta_eps, return_eps_ne=True
        )
        Q = (eps_rate_zeta - eps_rate_ref) / eps_zeta[:, None, None]

        # Ill-conditioning guard (e.g., at zero-stress far-field cells where
        # the residue is insensitive to ζ).  Zeroing the tangent contribution
        # for these cells is equivalent to making ζ locally inert for this
        # Newton iteration.
        H_MIN = 1e-12
        self.ind_h_small = to.where(to.abs(self.h) < H_MIN)[0]
        if len(self.ind_h_small) > 0:
            self.h[self.ind_h_small] = 1.0

        B = (self.r / self.h)[:, None, None] * Q

        # --- P = dr/dsigma via stress-perturbation ----------------------- #
        # Probe only the upper-triangular component of the symmetric tensor,
        # matching Desai: _compute_md_fields reads only (i,j) with i<=j, so
        # the FD is consistent with how the residue sees stress.  The factor
        # of 2 on shear entries is carried inside _compute_H via the Voigt
        # packing.
        self.P = to.zeros_like(stress)
        stress_eps = stress.clone()
        for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
            stress_eps[:, i, j] += EPSILON_STRESS
            r_sig = self.compute_residue(stress_eps, self.zeta, Temp, dt)
            self.P[:, i, j] = (r_sig - self.r) / EPSILON_STRESS
            self.P[:, j, i] = self.P[:, i, j]
            stress_eps[:, i, j] -= EPSILON_STRESS

        H = self._compute_H(Q, self.P)
        H_over_h = H / self.h[:, None, None]

        if len(self.ind_h_small) > 0:
            B[self.ind_h_small] = 0.0
            H_over_h[self.ind_h_small] = 0.0
            self.P[self.ind_h_small] = 0.0

        return B, H_over_h

    def _compute_H(self, Q: to.Tensor, P: to.Tensor) -> to.Tensor:
        """
        Tensorial-Voigt 6x6 packing of the rank-one tensor Q (outer) P.
        Shear rows/columns carry the factor of 2 implicit in tensorial
        Voigt storage of symmetric tensors — matches
        `ViscoplasticDesai.compute_H`.
        """
        n_elems = P.shape[0]
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
