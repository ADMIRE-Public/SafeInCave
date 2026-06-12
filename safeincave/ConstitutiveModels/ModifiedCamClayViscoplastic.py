# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
import torch as to
from .NonElasticElement import NonElasticElement


class ModifiedCamClayViscoplastic(NonElasticElement):
    """
    Modified Cam-Clay viscoplastic element with Perzyna overstress and
    theta-softening of the preconsolidation pressure p_c.

    Constitutive equations (formula sheet, compression-positive convention):
        p          = trace(sigma_pos) / 3
        s_ij       = sigma_pos_ij - p delta_ij
        q          = sqrt(3 J2) = sqrt(1.5 s:s)
        F          = p^2 - p p_c + q^2 / M^2
        lambda_P   = <F / p_c^2>^n / eta_v             (Macaulay bracket)
        epsdot_ij  = lambda_P * dF/dsigma_ij           (associated flow)
        p_n        = p + q^2 / (M^2 p)                 (image point on F)
        p_c_new    = p_c_old * (p_n / p_c_old)^theta   (cyclic-hardening)

    Numerical scheme
    ----------------
    p_c is treated as an internal state variable. The per-step update above
    is recast as a Newton residue and linearised by finite-difference
    probing of both p_c and sigma (same pattern as `MunsonDawsonCreep`):

        r(p_c, sigma) = p_c - p_c_old * (p_n(sigma) / p_c_old)^theta

    Parameters
    ----------
    M, lam, kap, theta, pc0, e0, eta_v, n_rate : torch.Tensor
        Per-element MCC parameters, shape (N,). ``lam > kap`` is required.
    name : str, optional
        Element identifier, default ``"cam_clay"``.

    Attributes
    ----------
    pc, pc_old : torch.Tensor, shape (N,)
        Current iterate and last-committed preconsolidation pressure [Pa].
    Fvp : torch.Tensor, shape (N,)
        Diagnostic: last evaluated yield-function value [Pa^2].
    r, h, P : torch.Tensor
        Residue and FD sensitivities populated by
        ``compute_B_and_H_over_h`` and consumed by
        ``increment_internal_variables``.

    Notes
    -----
    Adapted from MCC notebooks provided by Saeed
    (``examples/cam_clay/Calibration_Clay_Rich.ipynb``,
    ``Cyclic_MCCM_Clay_Rich.ipynb``) — parameter conventions and
    yield-surface form follow those scripts.
    """

    _SQRT_FLOAT64_EPS = 1.4901161193847656e-8  # math.sqrt(2.220446e-16)

    def __init__(
        self,
        M: to.Tensor,
        lam: to.Tensor,
        kap: to.Tensor,
        theta: to.Tensor,
        pc0: to.Tensor,
        e0: to.Tensor,
        eta_v: to.Tensor,
        n_rate: to.Tensor,
        name: str = "cam_clay",
    ):
        super().__init__(M.shape[0])
        self.name = name

        self.M = M.to(dtype=to.float64)
        self.lam = lam.to(dtype=to.float64)
        self.kap = kap.to(dtype=to.float64)
        self.theta = theta.to(dtype=to.float64)
        self.e0 = e0.to(dtype=to.float64)
        self.eta_v = eta_v.to(dtype=to.float64)
        self.n_rate = n_rate.to(dtype=to.float64)

        # Internal state: preconsolidation pressure
        self.pc = pc0.to(dtype=to.float64).clone()
        self.pc_old = self.pc.clone()

        # Diagnostic: last yield-function value
        self.Fvp = to.zeros(self.n_elems, dtype=to.float64)

        # Newton-coupling storage (filled by compute_B_and_H_over_h)
        self.r = to.zeros(self.n_elems, dtype=to.float64)
        self.h = to.ones(self.n_elems, dtype=to.float64)
        self.P = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.ind_h_small = to.tensor([], dtype=to.long)

    # ------------------------------------------------------------------ #
    # ISV lifecycle
    # ------------------------------------------------------------------ #

    def update_internal_variables(self) -> None:
        """Commit p_c at end of a converged time step."""
        self.pc_old = self.pc.clone()

    def increment_internal_variables(
        self, stress: to.Tensor, stress_k: to.Tensor, dt: float
    ) -> None:
        """
        Apply the linearised Newton correction for p_c using the stored
        residue and sensitivities from the most recent
        ``compute_B_and_H_over_h``:

            delta_pc = -(r + P : (stress - stress_k)) / h
        """
        delta_sigma = stress - stress_k
        delta_pc = -(self.r + to.einsum("bij,bij->b", self.P, delta_sigma)) / self.h
        if len(self.ind_h_small) > 0:
            delta_pc[self.ind_h_small] = 0.0
        self.pc = self.pc + delta_pc
        # p_c is physically positive — floor at 1 Pa as a numerical guard
        self.pc = to.clamp(self.pc, min=1.0)

    # ------------------------------------------------------------------ #
    # Physics evaluator (shared by rate, residue, and FD probes)
    # ------------------------------------------------------------------ #

    def _compute_mcc_fields(self, stress: to.Tensor, pc: to.Tensor):
        """
        Compute MCC intermediate quantities for a given (stress, p_c) state.

        Parameters
        ----------
        stress : (N, 3, 3) — safeincave convention (compression-negative, Pa)
        pc     : (N,)     — preconsolidation pressure (Pa, positive)

        Returns
        -------
        s_dev      : (N, 3, 3) deviatoric part of compression-positive stress
        p_safe     : (N,) mean effective stress, compression-positive, floored
        F          : (N,) yield-function value F = p^2 - p p_c + q^2 / M^2
        lambda_P   : (N,) Perzyna multiplier <F/p_c^2>^n / eta_v  [1/s]
        m_pos      : (N, 3, 3) dimensionless flow direction
                     m = (1/p_c) * dF/dsigma_ij in compression-positive frame
        p_n        : (N,) consistent preconsolidation pressure p + q^2/(M^2 p)
        """
        # Compression-positive stress
        s_pos = -stress
        s_xx = s_pos[:, 0, 0]
        s_yy = s_pos[:, 1, 1]
        s_zz = s_pos[:, 2, 2]
        s_xy = s_pos[:, 0, 1]
        s_xz = s_pos[:, 0, 2]
        s_yz = s_pos[:, 1, 2]

        p = (s_xx + s_yy + s_zz) / 3.0
        s_dev = s_pos.clone()
        s_dev[:, 0, 0] = s_xx - p
        s_dev[:, 1, 1] = s_yy - p
        s_dev[:, 2, 2] = s_zz - p

        # q^2 = 3 J2 = 1.5 s:s
        q2 = 1.5 * (
            s_dev[:, 0, 0] ** 2
            + s_dev[:, 1, 1] ** 2
            + s_dev[:, 2, 2] ** 2
            + 2.0 * (s_xy**2 + s_xz**2 + s_yz**2)
        )

        # Numerical floors: 1 Pa on p and pc to guard divisions / power laws
        p_safe = to.clamp(p, min=1.0)
        pc_safe = to.clamp(pc, min=1.0)

        M2 = self.M * self.M
        F = p_safe * p_safe - p_safe * pc_safe + q2 / M2

        # Perzyna multiplier: <F/pc^2>^n / eta_v
        eta_safe = to.clamp(self.eta_v, min=1e-30)
        ratio = to.clamp(F / (pc_safe * pc_safe), min=0.0)
        ratio = to.clamp(ratio, max=1e6)  # cap to avoid FD overflow
        lambda_P = (ratio**self.n_rate) / eta_safe

        # Dimensionless flow direction m = (1/pc) * dF/dsigma_pos.
        # dF/dsigma_pos = (2p - pc)/3 delta + (3/M^2) s_dev  [Pa]
        # so m_vol = (2p - pc) / (3 pc), m_dev = (3 / (M^2 pc)) s_dev (both dimensionless).
        coeff_vol = (2.0 * p_safe - pc_safe) / (3.0 * pc_safe)
        coeff_dev = 3.0 / (M2 * pc_safe)
        m_pos = coeff_dev[:, None, None] * s_dev
        m_pos[:, 0, 0] += coeff_vol
        m_pos[:, 1, 1] += coeff_vol
        m_pos[:, 2, 2] += coeff_vol

        # Consistent preconsolidation (root of F=0 in p_c): p_n = p + q^2/(M^2 p)
        p_n = p_safe + q2 / (M2 * p_safe)

        return s_dev, p_safe, F, lambda_P, m_pos, p_n

    def compute_residue(self, stress: to.Tensor, pc: to.Tensor, dt: float) -> to.Tensor:
        """
        Per-step residue for p_c using Saeed's cyclic-hardening update
        (Calibration_Clay_Rich.ipynb, ``simulate()``):

            p_c_new = p_c_old * (p_n / p_c_old) ** theta
            p_n     = p + q^2 / (M^2 p)

        Recast as a Newton residue:

            r(p_c, sigma) = p_c - p_c_old * (p_n(sigma) / p_c_old) ** theta

        ``dt`` is unused (Saeed's update is per stress increment, not per
        unit time) but kept in the signature for interface symmetry with
        :meth:`MunsonDawsonCreep.compute_residue`.
        """
        _, _, _, _, _, p_n = self._compute_mcc_fields(stress, pc)

        pc_old_safe = to.clamp(self.pc_old, min=1.0)
        p_n_safe = to.clamp(p_n, min=1.0)

        pc_target = pc_old_safe * (p_n_safe / pc_old_safe) ** self.theta
        return pc - pc_target

    # ------------------------------------------------------------------ #
    # Strain-rate evaluator
    # ------------------------------------------------------------------ #

    def compute_eps_ne_rate(
        self,
        stress: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        pc: to.Tensor = None,
        return_eps_ne: bool = False,
    ):
        """
        Viscoplastic strain rate at (sigma, p_c).

        Parameters
        ----------
        stress : (N, 3, 3) safeincave-convention stress (Pa)
        phi1   : unused (kept for interface compatibility)
        Temp   : unused (kept for interface compatibility)
        pc     : optional override for FD probing; defaults to ``self.pc``.
        return_eps_ne : if True, return rate without touching ``self.eps_ne_rate``
        """
        if pc is None:
            pc = self.pc

        _, _, F, lambda_P, m_pos, _ = self._compute_mcc_fields(stress, pc)

        # Convert compression-positive flow direction to safeincave convention.
        # m_pos is dimensionless, lambda_P is in 1/s, so eps_rate is in 1/s.
        eps_rate = -m_pos * lambda_P[:, None, None]

        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
            self.Fvp = F.clone()

    # ------------------------------------------------------------------ #
    # ISV-coupled consistent tangent
    # ------------------------------------------------------------------ #

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta_t: float, Temp: to.Tensor
    ):
        """
        FD-based ISV-coupled consistent-tangent terms (mirrors
        :meth:`MunsonDawsonCreep.compute_B_and_H_over_h`):

            r  = residue at (sigma, pc)
            h  = (r(sigma, pc+eps) - r) / eps
            Q  = (eps_vp(sigma, pc+eps) - eps_vp) / eps   (N, 3, 3)
            P  = (r(sigma+eps*e_ij, pc) - r) / eps        (N, 3, 3, symmetric)
            H  = Q (outer) P  packed in tensorial-Voigt 6x6
            B  = (r / h) * Q
        """
        # Forward-difference scale for pc: tied to pc itself + pc0 floor
        pc_scale = to.clamp(to.abs(self.pc), min=1.0)
        eps_pc = self._SQRT_FLOAT64_EPS * pc_scale  # (N,)

        # Forward-difference scale for stress (Pa) — same as MD
        EPSILON_STRESS = 1e-1

        # --- r, h, Q via pc perturbation -------------------------------- #
        self.r = self.compute_residue(stress, self.pc, dt)

        pc_eps = self.pc + eps_pc
        r_pc = self.compute_residue(stress, pc_eps, dt)
        self.h = (r_pc - self.r) / eps_pc

        eps_rate_ref = self.compute_eps_ne_rate(
            stress, dt * theta_t, Temp, pc=self.pc, return_eps_ne=True
        )
        eps_rate_pc = self.compute_eps_ne_rate(
            stress, dt * theta_t, Temp, pc=pc_eps, return_eps_ne=True
        )
        Q = (eps_rate_pc - eps_rate_ref) / eps_pc[:, None, None]

        H_MIN = 1e-12
        self.ind_h_small = to.where(to.abs(self.h) < H_MIN)[0]
        if len(self.ind_h_small) > 0:
            self.h[self.ind_h_small] = 1.0

        B = (self.r / self.h)[:, None, None] * Q

        # --- P = dr/dsigma via stress perturbation ---------------------- #
        self.P = to.zeros_like(stress)
        stress_eps = stress.clone()
        for i, j in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]:
            stress_eps[:, i, j] += EPSILON_STRESS
            r_sig = self.compute_residue(stress_eps, self.pc, dt)
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
        Tensorial-Voigt 6x6 packing of Q (outer) P. Identical to
        :meth:`MunsonDawsonCreep._compute_H` — shear rows/columns carry the
        factor of 2 for tensorial Voigt storage of symmetric tensors.
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
