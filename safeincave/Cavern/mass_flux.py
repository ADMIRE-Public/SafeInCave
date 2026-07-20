# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

import CoolProp.CoolProp as CP
import numpy as np
from ..Thermo.CavernThermodynamics import CavernThermodynamics
from .base import Cavern
from .geometry import CavernVolumeComputer, HeatFluxComputer


class Cavern_MassFlux(Cavern):
    """
    Cavern model advanced from a prescribed mass-flux history.

    The cavern mass is updated from the interpolated mass flux and time-step
    size. Pressure, temperature, and density are then solved with
    :class:`CavernThermodynamics`, including heat exchange and volume change.

    Parameters
    ----------
    grid : GridHandlerGMSH-like
        Grid object used for heat-flux and volume computations.
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    sym_scale : int
        Symmetry multiplier for volume computation. Valid values are ``1``,
        ``2``, and ``4``.
    reference_point : list of float, optional
        Point used by :class:`CavernVolumeComputer` for volume integration.
    fluid : str
        CoolProp fluid name.
    P_init : float
        Initial gauge pressure, in Pa.
    T_init : float
        Initial fluid temperature, in K.
    T_in : float
        Temperature of injected fluid, in K.
    Mflux_values : list of float
        Prescribed mass-flux values over time.
    time_values : list of float
        Times corresponding to ``Mflux_values``.
    direction : int, default=2
        Coordinate direction used to calculate the cavern midpoint.
    g : float, default=-9.81
        Gravitational acceleration.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.
    P_atm : float, default=101325.0
        Atmospheric pressure used to convert gauge pressure to absolute pressure.

    Attributes
    ----------
    type : str
        Always ``"Cavern_MFlux"``.
    Mflux : float
        Current mass-flux value.
    P : float
        Current gauge pressure.
    T : float
        Current temperature.
    density : float
        Current fluid density.
    Q : float
        Integrated heat over the current time step.
    heat : HeatFluxComputer
        Helper used to integrate heat flux on the cavern boundary.
    cvc : CavernVolumeComputer
        Helper used to compute cavern volume.
    model : CavernThermodynamics
        Thermodynamic model used to update pressure, temperature, and density.
    V_hist, M_hist, P_hist, T_hist, Q_hist, density_hist, t_hist : list
        Recorded histories for volume, mass, pressure, temperature, heat,
        density, and time.
    """

    def __init__(
        self,
        *,
        grid: any,
        cavern_name: str,
        sym_scale: int,
        reference_point: list[float, float, float] = None,
        fluid: str,
        P_init: float,  # Initial gauge pressure (Pa)
        T_init: float,  # Initial fluid temperature (K)
        T_in: float,  # Temperature of injected fluid (K)
        Mflux_values: list,
        time_values: list,
        direction: int = 2,
        g: float = -9.81,
        h_conv: float = None,
        P_atm: float = 101325.0,  # Atmospheric pressure in Pa
    ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_MFlux"
        self.Mflux_values = Mflux_values
        self.time_values = time_values
        self.Mflux = self.Mflux_values[0]
        self.direction = direction
        self.gravity = g
        self.P_atm = P_atm
        self.P_init = P_init
        self.T_init = T_init
        self.T_in = T_in
        self.Q = 0.0

        # Initial gauge pressure and temperature
        self.P = P_init
        self.T = T_init
        self.P0 = P_init
        self.T0 = T_init

        # Initialize heat flux computer
        self.heat = HeatFluxComputer(grid=grid, boundary_name=self.cavern_name)

        # Initialize cavern volume computer
        self.cvc = CavernVolumeComputer(
            grid=grid,
            boundary_name=self.cavern_name,
            reference_point=reference_point,
            sym_scale=sym_scale,
        )
        self.ref_pos = self.cvc.calculate_cavern_midpoint(direction=self.direction)

        # Initialize thermodynamic model
        self.model = CavernThermodynamics(self.fluid)

        # Calculate initial density
        AS = CP.AbstractState("HEOS", self.fluid)
        AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = AS.rhomass()

        # Initialize histories
        self.V_hist = []
        self.M_hist = []
        self.P_hist = []
        self.T_hist = []
        self.Q_hist = []
        self.density_hist = []
        self.t_hist = []

    def calculate_volume(self, u: any = None) -> None:
        if u is None:
            self.V = self.cvc.compute()
            self.ref_pos = self.cvc.calculate_cavern_midpoint(direction=self.direction)
        else:
            self.V = self.cvc.compute(u)
            self.ref_pos = self.cvc.calculate_cavern_midpoint(
                direction=self.direction, u=u
            )

    def calculate_heat(self, dt: float, T: any = None, kappa: any = None) -> None:
        self.Q = self.heat.compute(dt, T, kappa)
        # self.Q = 0.0

    def update_cavern(self, t: float, dt: float) -> None:
        Mflux = np.interp(t, self.time_values, self.Mflux_values)
        self.M = self.M0 + Mflux * dt
        dm = self.M - self.M0
        self.P, self.T, self.density = self.model.solve(
            dm=dm,
            Q_in=self.Q,
            T_in=self.T_in,
            P0=self.P0 + self.P_atm,  # Convert to absolute pressure
            T0=self.T0,
            V0=self.V0,
            V1=self.V,
        )
        if self.P <= 0.0:
            raise ValueError(f"P must be > 0, got {self.P}")
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")

        # Convert to gauge pressure
        self.P -= self.P_atm
        # print(self.P/1e6, self.T, self.density, self.M0, dm, self.V)

    def record_data(self, t: float) -> None:
        self.density_hist.append(self.density)
        self.V_hist.append(self.V)
        self.M_hist.append(self.V * self.density)
        self.P_hist.append(self.P)
        self.T_hist.append(self.T)
        self.Q_hist.append(self.Q)
        self.t_hist.append(t)

    def calculate_initial_condition(self):
        self.V0 = self.V
        self.P0 = self.P
        self.T0 = self.T
        self.M0 = self.V * self.density

    def create_data(self):
        data = {}
        data["Time"] = self.t_hist
        data["Pressure"] = self.P_hist
        data["Temperature"] = self.T_hist
        data["Density"] = self.density_hist
        data["Volume"] = self.V_hist
        data["Mass"] = self.M_hist
        data["Heat"] = self.Q_hist
        return data
