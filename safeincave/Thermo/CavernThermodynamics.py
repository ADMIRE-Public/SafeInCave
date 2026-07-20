# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from typing import Tuple
import CoolProp.CoolProp as CP
import numpy as np
from scipy.optimize import root


class CavernThermodynamics(object):
    """
    One-step cavern/tank update using Helmholtz EOS via CoolProp
    and scipy.optimize.root, with unknowns [log(P1), T1].

    Governing equations:
        F1 = M1*u(P1,T1) - M0*u(P0,T0) - Q + W - m_inj*h(Pin,Tin) + m_out*h(Pout,Tout) = 0
        F2 = rho(P1,T1) - M1/V1 = 0

    Sign convention:
        Q > 0 : heat added to the gas
        W > 0 : boundary work done by the gas
    """

    def __init__(self, fluid: str):
        self.fluid = fluid
        self.AS = CP.AbstractState("HEOS", fluid)

    def _update_PT(self, P: float, T: float) -> None:
        if P <= 0.0:
            raise ValueError(f"P must be > 0, got {P}")
        if T <= 0.0:
            raise ValueError(f"T must be > 0, got {T}")
        self.AS.update(CP.PT_INPUTS, P, T)

    def rho_u_h(self, P: float, T: float) -> Tuple[float, float, float]:
        self._update_PT(P, T)
        rho = self.AS.rhomass()
        u = self.AS.umass()
        h = self.AS.hmass()
        return rho, u, h

    def solve_withdrawal(
        self,
        *,
        P0: float,
        T0: float,
        V0: float,
        m_prod: float = 0.0,
        Q: float,
        V1: float,
        logP1_guess: float | None = None,
        T1_guess: float | None = None,
        method: str = "hybr",
        tol: float = 1e-8,
    ) -> Tuple[float, float, float]:
        """
        Solve one time step.

        Parameters
        ----------
        P0, T0, V0 : float
            Initial cavern pressure [Pa], temperature [K], volume [m^3]
        Pin, Tin : float
            Inlet pressure [Pa] and temperature [K]
        m_prod : float
            Mass withdrawn (produced) during the time step [kg]
        Q : float
            Heat added to gas during the step [J]
        W : float
            Boundary work done by gas during the step [J]
        V1 : float
            Final cavern volume [m^3]
        logP1_guess, T1_guess : optional
            Initial guesses for unknowns log(P1) and T1
        method : str
            Solver used by scipy.optimize.root
        tol : float
            Nonlinear solver tolerance
        """
        if V0 <= 0.0 or V1 <= 0.0:
            raise ValueError("V0 and V1 must be positive.")
        if m_prod > 0.0:
            raise ValueError("m_prod must be nonpositive.")

        # Initial cavern state
        rho0, u0, _ = self.rho_u_h(P0, T0)
        M0 = rho0 * V0
        M1 = M0 + m_prod

        # Initial guess
        if logP1_guess is None:
            P1_guess = max(1e3, P0 * (M1 / max(M0, 1e-30)) * (V0 / V1))
            logP1_guess = np.log(P1_guess)

        if T1_guess is None:
            T1_guess = max(50.0, T0)

        x0 = np.array([logP1_guess, T1_guess], dtype=float)

        rho_target = M1 / V1

        # Outlet enthalpy
        Pout_guess = (P0 + P1_guess) / 2
        Tout_guess = (T0 + T1_guess) / 2
        _, _, h_out = self.rho_u_h(Pout_guess, Tout_guess)

        # Initial guess for boundary work
        W = (V1 - V0) * (P0 + P1_guess) / 2

        # Characteristic scaling for residuals
        energy_scale = max(1.0, abs(M0 * u0) + abs(Q) + abs(W) + abs(m_prod * h_out))
        density_scale = max(1.0, abs(rho_target))

        def residuals(x: np.ndarray) -> np.ndarray:
            logP1, T1 = float(x[0]), float(x[1])
            P1 = np.exp(logP1)

            Tout = (T0 + T1) / 2
            Pout = (P0 + P1) / 2
            W = (V1 - V0) * (P0 + P1) / 2

            try:
                rho1, u1, _ = self.rho_u_h(P1, T1)
                _, _, h_out = self.rho_u_h(Pout, Tout)
            except Exception:
                # Return large residuals if CoolProp fails at this state
                return np.array([1e20, 1e20], dtype=float)

            F1 = M1 * u1 - M0 * u0 - Q + W - m_prod * h_out
            F2 = rho1 - rho_target

            return np.array([F1 / energy_scale, F2 / density_scale], dtype=float)

        sol = root(
            residuals,
            x0,
            method=method,
            tol=tol,
        )

        logP1, T1 = sol.x
        rho1 = M1 / V1
        return float(np.exp(logP1)), float(T1), float(rho1)

    def solve_injection(
        self,
        *,
        P0: float,
        T0: float,
        V0: float,
        Tin: float,
        m_inj: float,
        Q: float,
        V1: float,
        logP1_guess: float | None = None,
        T1_guess: float | None = None,
        method: str = "hybr",
        tol: float = 1e-8,
    ) -> Tuple[float, float, float]:
        """
        Solve one time step.

        Parameters
        ----------
        P0, T0, V0 : float
            Initial cavern pressure [Pa], temperature [K], volume [m^3]
        Tin : float
            Inlet temperature [K]
        m_inj : float
            Injected mass during the time step [kg]
        Q : float
            Heat added to gas during the step [J]
        V1 : float
            Final cavern volume [m^3]
        logP1_guess, T1_guess : optional
            Initial guesses for unknowns log(P1) and T1
        method : str
            Solver used by scipy.optimize.root
        tol : float
            Nonlinear solver tolerance
        """
        if V0 <= 0.0 or V1 <= 0.0:
            raise ValueError("V0 and V1 must be positive.")
        if P0 <= 0.0:
            raise ValueError("P0 and Pin must be positive.")
        if T0 <= 0.0 or Tin <= 0.0:
            raise ValueError("T0 and Tin must be positive.")
        if m_inj <= 0.0:
            raise ValueError("m_inj must be positive.")

        # Initial cavern state
        rho0, u0, _ = self.rho_u_h(P0, T0)
        M0 = rho0 * V0
        M1 = M0 + m_inj

        # Initial guess
        if logP1_guess is None:
            P1_guess = max(1e3, P0 * (M1 / max(M0, 1e-30)) * (V0 / V1))
            logP1_guess = np.log(P1_guess)

        if T1_guess is None:
            T1_guess = max(50.0, T0)

        Pin_guess = (P0 + P1_guess) / 2

        # Initial solution guess
        x0 = np.array([logP1_guess, T1_guess], dtype=float)

        # Target density
        rho_target = M1 / V1

        # Inlet enthalpy
        _, _, h_in = self.rho_u_h(Pin_guess, Tin)

        # Initial guess for boundary work
        W = Pin_guess * (V1 - V0)

        # Characteristic scaling for residuals
        energy_scale = max(1.0, abs(M0 * u0) + abs(Q) + abs(W) + abs(m_inj * h_in))
        density_scale = max(1.0, abs(rho_target))

        def residuals(x: np.ndarray) -> np.ndarray:
            logP1, T1 = float(x[0]), float(x[1])
            P1 = np.exp(logP1)

            Pin = (P0 + P1) / 2
            W = (V1 - V0) * (P0 + P1) / 2

            try:
                rho1, u1, _ = self.rho_u_h(P1, T1)
            except Exception:
                # Return large residuals if CoolProp fails at this state
                return np.array([1e20, 1e20], dtype=float)

            try:
                _, _, h_in = self.rho_u_h(Pin, Tin)
            except Exception:
                # Return large residuals if CoolProp fails at this state
                return np.array([1e20, 1e20], dtype=float)

            F1 = M1 * u1 - M0 * u0 - Q + W - m_inj * h_in
            F2 = rho1 - rho_target

            return np.array([F1 / energy_scale, F2 / density_scale], dtype=float)

        sol = root(
            residuals,
            x0,
            method=method,
            tol=tol,
        )

        if not sol.success:
            print("Solver failed:", sol.message)

        logP1, T1 = sol.x
        rho1 = M1 / V1
        return float(np.exp(logP1)), float(T1), float(rho1)

    def solve(
        self,
        dm: float,
        Q_in: float,
        T_in: float,
        P0: float,
        T0: float,
        V0: float,
        V1: float,
    ):
        # print(P0, T0, V0, dm)
        if dm > 0.0:
            P1, T1, rho1 = self.solve_injection(
                P0=P0,
                T0=T0,
                V0=V0,
                Tin=T_in,
                m_inj=dm,
                Q=Q_in,
                V1=V1,
            )
        else:
            P1, T1, rho1 = self.solve_withdrawal(
                P0=P0,
                T0=T0,
                V0=V0,
                m_prod=dm,
                Q=Q_in,
                V1=V1,
            )
        return P1, T1, rho1
