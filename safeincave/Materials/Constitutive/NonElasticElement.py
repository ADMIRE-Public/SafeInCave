# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any
import torch as to
from ...Utils import dotdot_torch
from ...Derivatives import (
    DerivativeEvaluator,
    resolve_derivative_evaluator,
    to_namespace,
)


class NonElasticElement(ABC):
    """
    Abstract base for inelastic mechanisms (e.g., viscoelasticity,
    dislocation creep, viscoplasticity). Provides common storage
    and utility updates for internal variables.

    Parameters
    ----------
    n_elems : int
        Number of elements.
    derivative_method : str or DerivativeEvaluator, optional
        Backend used to differentiate the strain rate with respect to stress:
        ``"finite_difference"`` (default), ``"torch_ad"`` or ``"jax_ad"``.
        `None` uses the global default — see
        :func:`safeincave.set_default_derivative_method`.

    Attributes
    ----------
    n_elems : int
        Number of elements.
    eps_ne_rate, eps_ne_rate_old : torch.Tensor
        Current and previous non-elastic strain rate, shape (N, 3, 3).
    eps_ne_old, eps_ne_k : torch.Tensor
        Non-elastic strain at old and current time, shape (N, 3, 3).
    B : torch.Tensor
        State variable term (N, 3, 3) assembled in `compute_G_B`.
    G : torch.Tensor
        Tangent-like operator (N, 6, 6) assembled in `compute_G_B`.
    derivative : DerivativeEvaluator
        Resolved derivative backend used by `compute_E` and, where relevant,
        by `compute_B_and_H_over_h`.
    """

    #: Finite-difference step (Pa) used by `compute_E` under the FD backend.
    EPSILON_STRESS_RATE: float = 1e-2

    def __init__(
        self,
        n_elems: int,
        derivative_method: "str | DerivativeEvaluator | None" = None,
    ) -> None:
        self.n_elems = n_elems
        self.derivative = resolve_derivative_evaluator(derivative_method)
        self.eps_ne_rate = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.eps_ne_rate_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.eps_ne_old = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.eps_ne_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        self.G = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

    @abstractmethod
    def compute_eps_ne_rate(
        self,
        stress_vec: to.Tensor,
        phi1: float,
        Temp: to.Tensor,
        return_eps_ne: bool = False,
    ) -> None:
        pass

    def increment_internal_variables(self, *args) -> None:
        pass

    def update_internal_variables(self, *args) -> None:
        pass

    def compute_eps_ne_k(self, phi1: float, phi2: float) -> None:
        """
        Predictor for non-elastic strain at the previous iteration k.

        Parameters
        ----------
        phi1 : float
            Typically `phi1=dt*theta`.
        phi2 : float
            Typically `phi2=dt*(1-theta)`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.eps_ne_k`.
        """
        self.eps_ne_k = (
            self.eps_ne_old + phi1 * self.eps_ne_rate_old + phi2 * self.eps_ne_rate
        )

    def update_eps_ne_old(
        self, stress: to.Tensor, stress_k: to.Tensor, phi2: float
    ) -> None:
        """
        Update non-elastic strain from previous time step.

        Parameters
        ----------
        stress : torch.Tensor
            Stress at current iteration k+1, shape (N, 3, 3).
        stress_k : torch.Tensor
            Stress at previous iteration k, shape (N, 3, 3).
        phi2 : float
            Typically `phi2=dt*(1-theta)`.

        Returns
        -------
        None

        Side Effects
        ------------
        Updates `self.eps_ne_old`.
        """
        self.eps_ne_old = (
            self.eps_ne_k
            + phi2 * dotdot_torch(self.G, stress - stress_k)
            - phi2 * self.B
        )

    def update_eps_ne_rate_old(self) -> None:
        """
        Update the current non-nelastic strain rate to the previous time step.

        Returns
        -------
        None
        """
        self.eps_ne_rate_old = self.eps_ne_rate.clone()

    def rate_fn(
        self, xp: Any, stress: Any, phi1: float, Temp: Any, **isv: Any
    ) -> Any:
        """
        Array-namespace-generic strain-rate kernel used by the derivative backends.

        Parameters
        ----------
        xp : module
            Array namespace: `torch` for the finite-difference and torch-AD
            backends, `jax.numpy` for the JAX-AD backend.
        stress : array
            Stress per element, shape (N, 3, 3).
        phi1 : float
            Time integration factor (`dt*theta`).
        Temp : array
            Temperature per element.
        **isv
            Optional internal-state-variable overrides (e.g. `alpha`, `zeta`,
            `pc`), forwarded to the concrete model.

        Returns
        -------
        array
            Strain rate, shape (N, 3, 3), in the namespace `xp`.

        Notes
        -----
        The default implementation simply forwards to `compute_eps_ne_rate`,
        which covers the finite-difference and torch-AD backends. Mechanisms
        that want the JAX backend must override this with a kernel written
        against `xp` only — see :mod:`safeincave.Derivatives.namespace`.
        """
        if xp is to:
            return self.compute_eps_ne_rate(
                stress, phi1, Temp, return_eps_ne=True, **isv
            )
        raise NotImplementedError(
            f"{type(self).__name__} does not provide an array-namespace-generic "
            "strain-rate kernel, so derivative_method='jax_ad' is unavailable "
            "for it. Use 'finite_difference' or 'torch_ad', or override "
            "`rate_fn`."
        )

    def compute_E(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> None:
        """
        The 6×6 operator E = d(eps_ne)/d(stress).

        Evaluated by `self.derivative`: central finite differences by default,
        or exact forward-mode AD when the mechanism was built with
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
            Temperature per element.

        Returns
        -------
        torch.Tensor
            Operator `E` with shape (N, 6, 6).
        """
        phi1 = dt * theta
        return self.derivative.d_tensor_d_stress(
            lambda xp, s: self.rate_fn(xp, s, phi1, self._cast(xp, Temp)),
            stress,
            step=self.EPSILON_STRESS_RATE,
        )

    @staticmethod
    def _cast(xp: Any, *tensors: Any) -> Any:
        """Move torch parameters into the namespace `xp` (identity for torch)."""
        return to_namespace(xp, *tensors)

    def compute_B_and_H_over_h(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> None:
        """
        Compute state variable term `B` and linearization term `H/h`.

        Parameters
        ----------
        stress : torch.Tensor
            Stress per element, shape (N, 3, 3).
        dt : float
            Time step.
        theta : float
            Time integration parameter.
        Temp : torch.Tensor
            Temperature per element.

        Returns
        -------
        B : torch.Tensor
            Driving term, shape (N, 3, 3).
        H_over_h : torch.Tensor
            Linearization ratio, shape (N, 6, 6).

        Notes
        -----
        Default implementation returns zeros; subclasses override.
        """
        B = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
        H_over_h = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
        return B, H_over_h

    def compute_G_B(
        self, stress: to.Tensor, dt: float, theta: float, Temp: to.Tensor
    ) -> None:
        """
        Assemble `G` and `B` for the element based on `E` and `H/h`.

        Parameters
        ----------
        stress : torch.Tensor
        dt : float
        theta : float
        Temp : torch.Tensor

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.B` and `self.G`.
        """
        self.B, H_over_h = self.compute_B_and_H_over_h(stress, dt, theta, Temp)
        E = self.compute_E(stress, dt, theta, Temp)
        self.G = E.clone() - H_over_h.clone()

    def compute_T_IT(self) -> None:
        """
        Build volumetric coupling tensors `T` (3×3) and `IT` (6×6) from `G`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.T` and `self.IT`.
        """
        self.T = to.zeros((self.n_elems, 3, 3))
        self.T[:, 0, 0] = self.G[:, 0, 0] + self.G[:, 1, 0] + self.G[:, 2, 0]
        self.T[:, 1, 1] = self.G[:, 0, 1] + self.G[:, 1, 1] + self.G[:, 2, 1]
        self.T[:, 2, 2] = self.G[:, 0, 2] + self.G[:, 1, 2] + self.G[:, 2, 2]
        self.T[:, 1, 0] = self.T[:, 0, 1] = (
            self.G[:, 0, 3] + self.G[:, 1, 3] + self.G[:, 2, 3]
        ) / 2
        self.T[:, 2, 0] = self.T[:, 0, 2] = (
            self.G[:, 0, 4] + self.G[:, 1, 4] + self.G[:, 2, 4]
        ) / 2
        self.T[:, 2, 1] = self.T[:, 1, 2] = (
            self.G[:, 0, 5] + self.G[:, 1, 5] + self.G[:, 2, 5]
        ) / 2

        self.IT = to.zeros((self.n_elems, 6, 6))
        self.IT[:, 0, 0] = self.T[:, 0, 0]
        self.IT[:, 0, 1] = self.T[:, 1, 1]
        self.IT[:, 0, 2] = self.T[:, 2, 2]
        self.IT[:, 0, 3] = self.T[:, 0, 1] + self.T[:, 1, 0]
        self.IT[:, 0, 4] = self.T[:, 0, 2] + self.T[:, 2, 0]
        self.IT[:, 0, 5] = self.T[:, 1, 2] + self.T[:, 2, 1]
        self.IT[:, 1, :] = self.IT[:, 0, :]
        self.IT[:, 2, :] = self.IT[:, 0, :]

    def compute_Bvol_Tvol(self) -> None:
        """
        Compute volumetric parts of `B` and `T`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.B_vol` and `self.T_vol` as traces of `B` and `T`.
        """
        self.T_vol = to.einsum("bii->b", self.T)
        self.B_vol = to.einsum("bii->b", self.B)

    def compute_Gtilde_Btilde(self) -> None:
        """
        Compute deviatoric parts `G_tilde` and `B_tilde`.

        Returns
        -------
        None

        Side Effects
        ------------
        Sets `self.G_tilde` and `self.B_tilde`.
        """
        I_3x3 = to.eye(3).expand(self.n_elems, -1, -1)
        self.G_tilde = self.G - self.IT / 3
        self.B_tilde = self.B - self.B_vol[:, None, None] * I_3x3 / 3
