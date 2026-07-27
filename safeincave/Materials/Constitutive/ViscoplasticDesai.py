# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from typing import Any
import torch as to
import numpy as np
from ...Utils import MPa
from ...Derivatives import build_sym3, pack_H_voigt, split_sym3
from .NonElasticElement import NonElasticElement

SQRT27 = float(np.sqrt(27))


def _invariants_kernel(xp: Any, s: tuple, sigma_t: Any) -> tuple:
    """
    Stress invariants from components already in the Desai convention
    (compression-positive, MPa).

    Returns `(s, I1, I2, I3, J2, J2_safe, J3, Sr, I1_star, J2_pos)`.

    `J2_safe` is `J2` clamped away from zero and `J2_pos` is the validity
    mask.  Every division by `J2` downstream must use `J2_safe` and be gated
    by `J2_pos`: dividing by the raw `J2` produces NaNs at zero-deviator
    points which, unlike the historical index-assignment masking, would
    poison the AD tangent even where the mask discards the value.
    """
    s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = s

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

    J2_pos = J2 > 0.0
    J2_safe = xp.where(J2_pos, J2, xp.ones_like(J2))
    Sr = xp.where(J2_pos, -(J3 * SQRT27) / (2 * J2_safe**1.5), xp.zeros_like(J2))

    I1_star = I1 + sigma_t
    return s, I1, I2, I3, J2, J2_safe, J3, Sr, I1_star, J2_pos


def _invariants_from_stress(xp: Any, stress: Any, sigma_t: Any) -> tuple:
    """`_invariants_kernel` from a stress tensor in the safeincave convention."""
    return _invariants_kernel(xp, split_sym3(-stress / MPa), sigma_t)


def _yield_kernel(
    xp: Any, alpha: Any, I1: Any, J2: Any, Sr: Any, p: dict
) -> Any:
    """Desai yield function `Fvp = J2 + F1 * F2**m`."""
    F1 = alpha * I1 ** p["n"] - p["gamma"] * I1**2
    F2 = xp.exp(p["beta_1"] * I1) - p["beta"] * Sr
    return J2 + F1 * F2 ** p["m"]


def _rate_kernel(xp: Any, stress: Any, alpha: Any, p: dict) -> tuple:
    """
    Namespace-generic Desai viscoplastic rate kernel.

    Returns `(eps_vp_rate, Fvp)`.  Pure and out-of-place, with the two
    non-smooth branches (`J2 <= 0` and the Macaulay bracket on `Fvp`)
    expressed as double-`where` guards so neither value nor derivative can
    see a NaN from the discarded branch.
    """
    s, I1, I2, I3, J2, J2_safe, J3, Sr, I1_star, J2_pos = _invariants_from_stress(
        xp, stress, p["sigma_t"]
    )
    s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = s

    Fvp = _yield_kernel(xp, alpha, I1_star, J2, Sr, p)

    # Flow direction, i.e. d(Fvp)/d(stress), derived analytically.
    F1 = -alpha * I1_star ** p["n"] + p["gamma"] * I1_star**2
    F2 = xp.exp(p["beta_1"] * I1_star) - p["beta"] * Sr
    dF1_dI1 = 2 * p["gamma"] * I1_star - p["n"] * alpha * I1_star ** (p["n"] - 1)
    dF2m_dI1 = (
        p["beta_1"] * p["m"] * xp.exp(p["beta_1"] * I1_star) * F2 ** (p["m"] - 1)
    )
    dF_dI1 = -(dF1_dI1 * F2 ** p["m"] + F1 * dF2m_dI1)

    dF2_dJ2 = -(3 * p["beta"] * J3 * SQRT27) / (4 * J2_safe ** (5 / 2))
    dF_dJ2 = 1 - F1 * p["m"] * F2 ** (p["m"] - 1) * dF2_dJ2
    dF_dJ3 = (
        -p["m"] * F1 * p["beta"] * SQRT27 * F2 ** (p["m"] - 1) / (2 * J2_safe**1.5)
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

    dQdS = build_sym3(xp, dQdS_00, dQdS_11, dQdS_22, dQdS_01, dQdS_02, dQdS_12)
    # Wherever J2 <= 0, viscoplasticity is switched off.
    dQdS = xp.where(J2_pos[:, None, None], dQdS, xp.zeros_like(dQdS))

    # Macaulay bracket on Fvp, guarded so the inactive branch never evaluates
    # a fractional power of a negative number.
    yielding = Fvp > 0.0
    Fvp_safe = xp.where(yielding, Fvp, xp.ones_like(Fvp))
    lmbda = xp.where(
        yielding,
        p["mu_1"] * (Fvp_safe / p["F_0"]) ** p["N_1"],
        xp.zeros_like(Fvp),
    )

    return -dQdS * lmbda[:, None, None], Fvp


def _residue_kernel(
    xp: Any, eps_rate: Any, alpha: Any, qsi_old: Any, dt: float, p: dict
) -> tuple:
    """
    Hardening residue and the accumulated measure it depends on.

    Returns `(qsi, r)` with
    ``r = alpha - a_1 / (((a_1/alpha_0)**(1/eta) + qsi)**eta)``.

    The Frobenius norm is guarded: elastic elements have an exactly zero
    strain rate, where ``d(sqrt(u))/dx = u'/(2 sqrt(u))`` evaluates to 0/0.
    Finite differences never see this (they difference two zeros), but AD
    does, so the zero branch is short-circuited to a zero derivative — the
    same value the difference quotient produces.
    """
    sq = xp.sum(eps_rate**2, (-2, -1))
    positive = sq > 0.0
    sq_safe = xp.where(positive, sq, xp.ones_like(sq))
    norm = xp.where(positive, sq_safe**0.5, xp.zeros_like(sq))
    qsi = qsi_old + norm * dt
    r = alpha - p["a_1"] / (
        ((p["a_1"] / p["alpha_0"]) ** (1 / p["eta"]) + qsi) ** p["eta"]
    )
    return qsi, r


class ViscoplasticDesai(NonElasticElement):
    """
    Viscoplastic model of Desai-type with hardening (state) variable `alpha`.

    Parameters
    ----------
    mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0 : torch.Tensor
        Model parameters per element, shape (N,).
    name : str, optional
        Element name, by default "desai".
    derivative_method : str or DerivativeEvaluator, optional
        Derivative backend used for both `compute_E` and the ISV-coupled
        terms in `compute_B_and_H_over_h`; see :class:`NonElasticElement`.

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
        derivative_method: Any = None,
    ) -> None:
        super().__init__(mu_1.shape[0], derivative_method=derivative_method)
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

    _PARAM_NAMES = (
        "mu_1",
        "N_1",
        "a_1",
        "eta",
        "n",
        "beta_1",
        "beta",
        "m",
        "gamma",
        "sigma_t",
        "alpha_0",
    )

    def _params(self, xp: Any) -> dict:
        """Material parameters cast into the array namespace `xp`."""
        values = self._cast(xp, *(getattr(self, k) for k in self._PARAM_NAMES))
        params = dict(zip(self._PARAM_NAMES, values))
        params["F_0"] = self.F_0
        return params

    def rate_fn(
        self, xp: Any, stress: Any, phi1: float, Temp: Any, alpha: Any = None
    ) -> Any:
        """Namespace-generic strain-rate kernel (see :class:`NonElasticElement`)."""
        if alpha is None:
            alpha = self._cast(xp, self.alpha)
        eps_rate, _ = _rate_kernel(xp, stress, alpha, self._params(xp))
        return eps_rate

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
        self.qsi, r = _residue_kernel(
            to, eps_rate, alpha, self.qsi_old, dt, self._params(to)
        )
        return r

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
        # The components already are in the Desai convention
        # (compression-positive, MPa), so they go straight into the kernel.
        _, I1, I2, I3, J2, _, J3, Sr, I1_star, J2_pos = _invariants_kernel(
            to, (s_xx, s_yy, s_zz, s_xy, s_xz, s_yz), self.sigma_t
        )
        return I1, I2, I3, J2, J3, Sr, I1_star, to.where(~J2_pos)[0]

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

    def compute_Fvp(self, alpha: to.Tensor, I1: to.Tensor, J2: to.Tensor, Sr: to.Tensor) -> to.Tensor:
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
        return _yield_kernel(to, alpha, I1, J2, Sr, self._params(to))

    def compute_initial_hardening(self, stress: to.Tensor, Fvp_0: float = 0.0) -> None:
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
        alpha: to.Tensor | None = None,
        return_eps_ne: bool = False,
    ) -> None:
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

        eps_vp_rate, Fvp = _rate_kernel(to, stress, alpha, self._params(to))

        if return_eps_ne:
            return eps_vp_rate
        self.Fvp = Fvp.clone()
        self.eps_ne_rate = eps_vp_rate

    #: Forward-difference step (Pa) for dr/dsigma under the FD backend.
    EPSILON_STRESS_RESIDUE: float = 1e-1
    #: Relative forward-difference step for d/dalpha under the FD backend.
    EPSILON_ALPHA_REL: float = 1e-4

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> tuple[to.Tensor, to.Tensor]:
        """
        Compute `B` and `H/h` from the sensitivities of the hardening residue.

            r  = residue at the reference (sigma, alpha)
            h  = dr/dalpha                                    (N,)
            Q  = d(eps_vp)/dalpha                             (N, 3, 3)
            P  = dr/dsigma                                    (N, 3, 3, symmetric)
            H  = Q (outer) P  packed in tensorial Voigt 6x6
            B  = (r / h) * Q

        The three derivatives come from `self.derivative`: forward finite
        differences by default, exact AD when the element was built with
        ``derivative_method="torch_ad"`` / ``"jax_ad"``.

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
        The residue `r` — the Newton linearisation point — is referenced to
        the last computed strain rate `self.eps_ne_rate`, while the
        accumulated measure `self.qsi` is evaluated from the rate at the
        *current* stress.  Both match the previous behaviour: `qsi` used to
        be a side effect of `compute_residue` and ended up holding whatever
        the last finite-difference probe left behind, which was the rate at
        a stress perturbed by 1e-1 Pa.  Evaluating it at the unperturbed
        stress gives the same number without depending on the differentiation
        scheme.
        """
        params = self._params(to)
        eps_alpha = self.EPSILON_ALPHA_REL * self.alpha

        # Newton linearisation point: residue at the last computed rate.
        _, self.r = _residue_kernel(
            to, self.eps_ne_rate, self.alpha, self.qsi_old, dt, params
        )

        # Accumulated inelastic measure, from the rate at the current stress.
        rate_now, _ = _rate_kernel(to, stress, self.alpha, params)
        self.qsi, _ = _residue_kernel(
            to, rate_now, self.alpha, self.qsi_old, dt, params
        )

        def _isv_probe(xp, alpha):
            p = self._params(xp)
            rate, _ = _rate_kernel(xp, self._cast(xp, stress), alpha, p)
            _, r = _residue_kernel(
                xp, rate, alpha, self._cast(xp, self.qsi_old), dt, p
            )
            return rate, r

        Q, self.h = self.derivative.d_d_isv(
            _isv_probe,
            self.alpha,
            step=eps_alpha,
            f0=(self.eps_ne_rate, self.r),
        )

        B = (self.r / self.h)[:, None, None] * Q

        def _stress_probe(xp, s):
            p = self._params(xp)
            alpha = self._cast(xp, self.alpha)
            rate, _ = _rate_kernel(xp, s, alpha, p)
            _, r = _residue_kernel(
                xp, rate, alpha, self._cast(xp, self.qsi_old), dt, p
            )
            return r

        self.P = self.derivative.d_scalar_d_stress(
            _stress_probe,
            stress,
            step=self.EPSILON_STRESS_RESIDUE,
            f0=self.r,
        )

        H_over_h = pack_H_voigt(Q, self.P) / self.h[:, None, None]

        return B, H_over_h
