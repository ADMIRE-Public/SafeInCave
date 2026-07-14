# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

import CoolProp.CoolProp as CP
import numpy as np
from .base import Cavern
from .geometry import CavernVolumeComputer, HeatFluxComputer


class Cavern_PT(Cavern):
    """
    Cavern model with prescribed time-dependent pressure and temperature.

    Pressure and temperature are interpolated from user-provided histories.
    The current density is evaluated with CoolProp using the absolute pressure
    ``P + P_atm`` and temperature ``T``.

    Parameters
    ----------
    grid : GridHandlerGMSH-like
        Grid object used for heat-flux and volume computations.
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    fluid : str
        CoolProp fluid name.
    sym_scale : int
        Symmetry multiplier for volume computation. Valid values are ``1``,
        ``2``, and ``4``.
    reference_point : list of float, optional
        Point used by :class:`CavernVolumeComputer` for volume integration.
    P_values : list of float
        Prescribed gauge pressure values over time, in Pa.
    T_values : list of float
        Prescribed temperature values over time, in K.
    time_values : list of float
        Times corresponding to ``P_values`` and ``T_values``.
    ref_pos : float, default=0.0
        Reference coordinate used by coupled cavern models.
    direction : int, default=0
        Coordinate direction associated with ``ref_pos`` and gravity.
    g : float, default=-9.81
        Gravitational acceleration.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.
    P_atm : float, default=101325.0
        Atmospheric pressure added to gauge pressure before CoolProp calls.

    Attributes
    ----------
    type : str
        Always ``"Cavern_PT"``.
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
    V_hist, M_hist, P_hist, T_hist, Q_hist, density_hist, t_hist : list
        Recorded histories for volume, mass, pressure, temperature, heat,
        density, and time.
    """

    def __init__(
        self,
        *,
        grid: any,
        cavern_name: str,
        fluid: str,
        sym_scale: int,
        reference_point: list[float, float, float] = None,
        P_values: list,  # Gauge pressure values (Pa)
        T_values: list,  # Temperature values (K)
        time_values: list,
        ref_pos: float = 0.0,
        direction: int = 0,
        g: float = -9.81,
        h_conv: float = None,
        P_atm: float = 101325.0,  # Atmospheric pressure in Pa
    ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_PT"
        self.P_values = P_values
        self.T_values = T_values
        self.time_values = time_values
        self.P_atm = P_atm
        self.ref_pos = ref_pos
        self.direction = direction
        self.gravity = g
        self.Q = 0.0

        self.AS = CP.AbstractState("HEOS", self.fluid)
        self.P = self.P_values[0]
        self.T = self.T_values[0]
        self.AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = self.AS.rhomass()

        # Initialize heat flux computer
        self.heat = HeatFluxComputer(grid=grid, boundary_name=self.cavern_name)

        # Initialize cavern volume computer
        self.cvc = CavernVolumeComputer(
            grid=grid,
            boundary_name=self.cavern_name,
            reference_point=reference_point,
            sym_scale=sym_scale,
        )

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
        else:
            self.V = self.cvc.compute(u)

    def calculate_heat(self, dt: float, T: any = None, kappa: any = None) -> None:
        self.Q = self.heat.compute(dt, T, kappa)

    def update_cavern(self, t: float, dt: float) -> None:
        self.P = np.interp(t, self.time_values, self.P_values)
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.P + self.P_atm <= 0.0:
            raise ValueError(
                f"Absolute pressure must be > 0, got {self.P + self.P_atm}"
            )
        if self.T <= 0.0:
            raise ValueError(f"Temperature must be > 0, got {self.T}")
        self.AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = self.AS.rhomass()

    def record_data(self, t: float) -> None:
        self.density_hist.append(self.density)
        self.V_hist.append(self.V)
        self.M_hist.append(self.density * self.V)
        self.P_hist.append(self.P)
        self.T_hist.append(self.T)
        self.Q_hist.append(self.Q)
        self.t_hist.append(t)

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
