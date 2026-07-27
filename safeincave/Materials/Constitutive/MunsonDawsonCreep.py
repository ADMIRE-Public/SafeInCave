# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
from ...Derivatives import deviator, pack_H_voigt, split_sym3
from .NonElasticElement import NonElasticElement


def _md_fields_kernel(
    xp: Any, stress: Any, Temp: Any, zeta: Any, p: dict
) -> tuple:
    """
    Namespace-generic Munson-Dawson intermediate quantities.

    Returns `(s_dev, sigma_safe, epsdot_ss, eps_t_star, F)`.  Pure and
    out-of-place: the piecewise hardening/recovery split is a `where` on
    both branches rather than boolean-mask assignment, so it traces cleanly
    and neither branch can leak a NaN into the derivative.
    """
    s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = split_sym3(stress)

    s_dev = deviator(xp, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)

    # sigma = sqrt(3 J2) = von Mises equivalent stress (Pa)
    sigma = xp.sqrt(
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
    sigma_safe = xp.clip(sigma, 1.0, None)
    mu_safe = xp.clip(p["mu"], 1.0, None)

    epsdot_ss = p["A"] * xp.exp(-p["Q"] / (p["R"] * Temp)) * (sigma_safe ** p["n"])

    ratio = xp.clip(sigma_safe / mu_safe, 1e-30, None)
    eps_t_star = p["K0"] * xp.exp(p["c"] * Temp) * (ratio ** p["m"])
    eps_t_star = xp.clip(eps_t_star, 1e-50, None)

    Delta_cap = p["alpha_w"] + p["beta_w"] * xp.log10(ratio)

    # Piecewise F with hardening / recovery branches.  The ±50 exponent
    # clamp guards against overflow at pathological stress states
    # (exp(50) ≈ 5e21 — far beyond any realistic F).
    r_arg = 1.0 - (zeta / eps_t_star)
    r_arg2 = r_arg * r_arg

    hardening = zeta <= eps_t_star
    exp_arg = xp.where(hardening, Delta_cap * r_arg2, -p["delta"] * r_arg2)
    F = xp.exp(xp.clip(exp_arg, -50.0, 50.0))

    return s_dev, sigma_safe, epsdot_ss, eps_t_star, F


def _rate_from_fields(s_dev: Any, sigma_safe: Any, epsdot_ss: Any, F: Any) -> Any:
    """Assemble ``epsdot_ij = F * epsdot_ss * (3/2) s_ij / sigma`` from fields."""
    scalar_rate = F * epsdot_ss
    flow_dir = (1.5 / sigma_safe)[:, None, None] * s_dev
    return flow_dir * scalar_rate[:, None, None]


def _rate_kernel(xp: Any, stress: Any, Temp: Any, zeta: Any, p: dict) -> Any:
    """Namespace-generic Munson-Dawson strain-rate kernel."""
    s_dev, sigma_safe, epsdot_ss, _, F = _md_fields_kernel(xp, stress, Temp, zeta, p)
    return _rate_from_fields(s_dev, sigma_safe, epsdot_ss, F)


def _residue_kernel(
    xp: Any, stress: Any, Temp: Any, zeta: Any, zeta_old: Any, dt: float, p: dict
) -> Any:
    """Backward-Euler residue ``r = zeta - zeta_old - (F - 1) * epsdot_ss * dt``."""
    _, _, epsdot_ss, _, F = _md_fields_kernel(xp, stress, Temp, zeta, p)
    return zeta - zeta_old - (F - 1.0) * epsdot_ss * dt


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
        derivative_method: Any = None,
    ) -> None:
        super().__init__(A.shape[0], derivative_method=derivative_method)
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

    _PARAM_NAMES = (
        "A",
        "Q",
        "n",
        "K0",
        "c",
        "m",
        "alpha_w",
        "beta_w",
        "delta",
        "mu",
    )

    def _params(self, xp: Any) -> dict:
        """Material parameters cast into the array namespace `xp`."""
        values = self._cast(xp, *(getattr(self, k) for k in self._PARAM_NAMES))
        params = dict(zip(self._PARAM_NAMES, values))
        params["R"] = self.R
        return params

    def rate_fn(
        self, xp: Any, stress: Any, phi1: float, Temp: Any, zeta: Any = None
    ) -> Any:
        """Namespace-generic strain-rate kernel (see :class:`NonElasticElement`)."""
        if zeta is None:
            zeta = self._cast(xp, self.zeta)
        return _rate_kernel(xp, stress, Temp, zeta, self._params(xp))

    def _compute_md_fields(
        self, stress_vec: to.Tensor, Temp: to.Tensor, zeta: to.Tensor
    ) -> tuple:
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
        return _md_fields_kernel(to, stress_vec, Temp, zeta, self._params(to))

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
        return _residue_kernel(
            to, stress, Temp, zeta, self.zeta_old, dt, self._params(to)
        )

    # ------------------------------------------------------------------ #
    # Strain-rate evaluator
    # ------------------------------------------------------------------ #

    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        zeta: to.Tensor | None = None,
        return_eps_ne: bool = False,
    ) -> None:
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

        s_dev, sigma_safe, epsdot_ss, eps_t_star, F = _md_fields_kernel(
            to, stress_vec, Temp, zeta, self._params(to)
        )
        eps_rate = _rate_from_fields(s_dev, sigma_safe, epsdot_ss, F)

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

    #: Forward-difference step (Pa) for dr/dsigma under the FD backend.
    EPSILON_STRESS_RESIDUE: float = 1e-1
    #: Elements whose |dr/dzeta| falls below this are treated as zeta-inert.
    H_MIN: float = 1e-12

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> tuple[to.Tensor, to.Tensor]:
        """
        Build the ISV-coupled consistent-tangent terms.

            r  = residue at the reference (sigma, zeta)
            h  = dr/dzeta                                     (N,)
            Q  = d(eps_MD)/dzeta                              (N, 3, 3)
            P  = dr/dsigma                                    (N, 3, 3, symmetric)
            H  = Q (outer) P  packed in tensorial-Voigt 6x6
            B  = (r / h) * Q

        The derivatives come from `self.derivative`: forward finite
        differences by default, exact AD when the element was built with
        ``derivative_method="torch_ad"`` / ``"jax_ad"``.

        r, h, P are stored on the instance for the subsequent
        `increment_internal_variables` call.  B and H/h are returned.
        """
        params = self._params(to)

        # Forward-difference scale for zeta: tied to the natural range eps_t*
        # so it stays meaningful whether zeta is 0 (start of run) or near eps_t*.
        _, _, _, eps_t_star_now, _ = _md_fields_kernel(
            to, stress, Temp, self.zeta, params
        )
        zeta_scale = to.clamp(to.abs(self.zeta) + eps_t_star_now, min=1e-30)
        eps_zeta = self._SQRT_FLOAT64_EPS * zeta_scale  # shape (N,)

        def _isv_probe(xp, zeta):
            p = self._params(xp)
            s = self._cast(xp, stress)
            T = self._cast(xp, Temp)
            rate = _rate_kernel(xp, s, T, zeta, p)
            r = _residue_kernel(
                xp, s, T, zeta, self._cast(xp, self.zeta_old), dt, p
            )
            return rate, r

        Q, self.h = self.derivative.d_d_isv(_isv_probe, self.zeta, step=eps_zeta)
        self.r = _residue_kernel(
            to, stress, Temp, self.zeta, self.zeta_old, dt, params
        )

        # Ill-conditioning guard (e.g., at zero-stress far-field cells where
        # the residue is insensitive to zeta).  Zeroing the tangent contribution
        # for these cells is equivalent to making zeta locally inert for this
        # Newton iteration.
        self.ind_h_small = to.where(to.abs(self.h) < self.H_MIN)[0]
        if len(self.ind_h_small) > 0:
            self.h[self.ind_h_small] = 1.0

        B = (self.r / self.h)[:, None, None] * Q

        # P = dr/dsigma.  Only the upper triangle is probed and then mirrored;
        # the factor of 2 on shear entries is carried by `pack_H_voigt`.
        def _stress_probe(xp, s):
            return _residue_kernel(
                xp,
                s,
                self._cast(xp, Temp),
                self._cast(xp, self.zeta),
                self._cast(xp, self.zeta_old),
                dt,
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
