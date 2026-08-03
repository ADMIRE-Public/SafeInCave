# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
from ...Derivatives import build_sym3, pack_H_voigt, split_sym3
from .NonElasticElement import NonElasticElement


def _mcc_fields_kernel(xp: Any, stress: Any, pc: Any, p: dict) -> tuple:
    """
    Namespace-generic Modified Cam-Clay intermediate quantities.

    Returns `(s_dev, p_safe, F, lambda_P, m_pos, p_n)`.  Pure and
    out-of-place; the Macaulay bracket on `F/pc^2` is a double-`where` so the
    inactive branch never evaluates a fractional power of zero, whose
    derivative would be infinite.
    """
    # Compression-positive stress
    s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = split_sym3(-stress)

    p_mean = (s_xx + s_yy + s_zz) / 3.0
    d_xx = s_xx - p_mean
    d_yy = s_yy - p_mean
    d_zz = s_zz - p_mean
    s_dev = build_sym3(xp, d_xx, d_yy, d_zz, s_xy, s_xz, s_yz)

    # q^2 = 3 J2 = 1.5 s:s
    q2 = 1.5 * (
        d_xx**2 + d_yy**2 + d_zz**2 + 2.0 * (s_xy**2 + s_xz**2 + s_yz**2)
    )

    # Numerical floors: 1 Pa on p and pc to guard divisions / power laws
    p_safe = xp.clip(p_mean, 1.0, None)
    pc_safe = xp.clip(pc, 1.0, None)

    M2 = p["M"] * p["M"]
    F = p_safe * p_safe - p_safe * pc_safe + q2 / M2

    # Perzyna multiplier: <F/pc^2>^n / eta_v
    eta_safe = xp.clip(p["eta_v"], 1e-30, None)
    ratio = xp.clip(F / (pc_safe * pc_safe), 0.0, 1e6)  # cap avoids overflow
    active = ratio > 0.0
    ratio_safe = xp.where(active, ratio, xp.ones_like(ratio))
    lambda_P = xp.where(
        active, (ratio_safe ** p["n_rate"]) / eta_safe, xp.zeros_like(ratio)
    )

    # Dimensionless flow direction m = (1/pc) * dF/dsigma_pos.
    # dF/dsigma_pos = (2p - pc)/3 delta + (3/M^2) s_dev  [Pa]
    # so m_vol = (2p - pc) / (3 pc), m_dev = (3 / (M^2 pc)) s_dev.
    coeff_vol = (2.0 * p_safe - pc_safe) / (3.0 * pc_safe)
    coeff_dev = 3.0 / (M2 * pc_safe)
    m_pos = build_sym3(
        xp,
        coeff_dev * d_xx + coeff_vol,
        coeff_dev * d_yy + coeff_vol,
        coeff_dev * d_zz + coeff_vol,
        coeff_dev * s_xy,
        coeff_dev * s_xz,
        coeff_dev * s_yz,
    )

    # Consistent preconsolidation (root of F=0 in p_c): p_n = p + q^2/(M^2 p)
    p_n = p_safe + q2 / (M2 * p_safe)

    return s_dev, p_safe, F, lambda_P, m_pos, p_n


def _rate_kernel(xp: Any, stress: Any, pc: Any, p: dict) -> tuple:
    """Namespace-generic MCC viscoplastic rate kernel; returns `(eps_rate, F)`."""
    _, _, F, lambda_P, m_pos, _ = _mcc_fields_kernel(xp, stress, pc, p)
    # m_pos is dimensionless and lambda_P is in 1/s, so eps_rate is in 1/s.
    # The sign converts the compression-positive direction to the
    # safeincave convention.
    return -m_pos * lambda_P[:, None, None], F


def _residue_kernel(
    xp: Any, stress: Any, pc: Any, pc_old: Any, p: dict
) -> Any:
    """Residue ``r = pc - pc_old * (p_n / pc_old)**theta``."""
    _, _, _, _, _, p_n = _mcc_fields_kernel(xp, stress, pc, p)

    pc_old_safe = xp.clip(pc_old, 1.0, None)
    p_n_safe = xp.clip(p_n, 1.0, None)

    pc_target = pc_old_safe * (p_n_safe / pc_old_safe) ** p["theta"]
    return pc - pc_target


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
        derivative_method: Any = None,
    ) -> None:
        super().__init__(M.shape[0], derivative_method=derivative_method)
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

    _PARAM_NAMES = ("M", "lam", "kap", "theta", "e0", "eta_v", "n_rate")

    def _params(self, xp: Any) -> dict:
        """Material parameters cast into the array namespace `xp`."""
        values = self._cast(xp, *(getattr(self, k) for k in self._PARAM_NAMES))
        return dict(zip(self._PARAM_NAMES, values))

    def rate_fn(
        self, xp: Any, stress: Any, phi1: float, Temp: Any, pc: Any = None
    ) -> Any:
        """Namespace-generic strain-rate kernel (see :class:`NonElasticElement`)."""
        if pc is None:
            pc = self._cast(xp, self.pc)
        eps_rate, _ = _rate_kernel(xp, stress, pc, self._params(xp))
        return eps_rate

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
        return _mcc_fields_kernel(to, stress, pc, self._params(to))

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
        return _residue_kernel(
            to, stress, pc, self.pc_old, self._params(to)
        )

    # ------------------------------------------------------------------ #
    # Strain-rate evaluator
    # ------------------------------------------------------------------ #

    def compute_eps_ne_rate(
        self,
        stress: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        pc: to.Tensor | None = None,
        return_eps_ne: bool = False,
    ) -> None:
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

        eps_rate, F = _rate_kernel(to, stress, pc, self._params(to))

        if return_eps_ne:
            return eps_rate
        else:
            self.eps_ne_rate = eps_rate
            self.Fvp = F.clone()

    # ------------------------------------------------------------------ #
    # ISV-coupled consistent tangent
    # ------------------------------------------------------------------ #

    #: Forward-difference step (Pa) for dr/dsigma under the FD backend.
    EPSILON_STRESS_RESIDUE: float = 1e-1
    #: Elements whose |dr/dpc| falls below this are treated as pc-inert.
    H_MIN: float = 1e-12

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta_t: float, Temp: to.Tensor
    ) -> tuple[to.Tensor, to.Tensor]:
        """
        ISV-coupled consistent-tangent terms (mirrors
        :meth:`MunsonDawsonCreep.compute_B_and_H_over_h`):

            r  = residue at (sigma, pc)
            h  = dr/dpc                                  (N,)
            Q  = d(eps_vp)/dpc                           (N, 3, 3)
            P  = dr/dsigma                               (N, 3, 3, symmetric)
            H  = Q (outer) P  packed in tensorial-Voigt 6x6
            B  = (r / h) * Q

        The derivatives come from `self.derivative`: forward finite
        differences by default, exact AD when the element was built with
        ``derivative_method="torch_ad"`` / ``"jax_ad"``.
        """
        params = self._params(to)

        # Forward-difference scale for pc: tied to pc itself + pc0 floor
        pc_scale = to.clamp(to.abs(self.pc), min=1.0)
        eps_pc = self._SQRT_FLOAT64_EPS * pc_scale  # (N,)

        def _isv_probe(xp, pc):
            p = self._params(xp)
            s = self._cast(xp, stress)
            rate, _ = _rate_kernel(xp, s, pc, p)
            r = _residue_kernel(xp, s, pc, self._cast(xp, self.pc_old), p)
            return rate, r

        Q, self.h = self.derivative.d_d_isv(_isv_probe, self.pc, step=eps_pc)
        self.r = _residue_kernel(to, stress, self.pc, self.pc_old, params)

        self.ind_h_small = to.where(to.abs(self.h) < self.H_MIN)[0]
        if len(self.ind_h_small) > 0:
            self.h[self.ind_h_small] = 1.0

        B = (self.r / self.h)[:, None, None] * Q

        def _stress_probe(xp, s):
            return _residue_kernel(
                xp,
                s,
                self._cast(xp, self.pc),
                self._cast(xp, self.pc_old),
                self._params(xp),
            )

        self.P = self.derivative.d_scalar_d_stress(
            _stress_probe,
            stress,
            step=self.EPSILON_STRESS_RESIDUE,
            f0=self.r,
        )

        H_over_h = pack_H_voigt(Q, self.P) / self.h[:, None, None]

        if len(self.ind_h_small) > 0:
            B[self.ind_h_small] = 0.0
            H_over_h[self.ind_h_small] = 0.0
            self.P[self.ind_h_small] = 0.0

        return B, H_over_h
