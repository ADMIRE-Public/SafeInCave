# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from .base import Cavern


class Cavern_T(Cavern):
    """
    Cavern model with prescribed time-dependent temperature.

    Parameters
    ----------
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    T_values : list of float
        Prescribed cavern temperatures over time.
    time_values : list of float
        Times corresponding to ``T_values``.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.

    Attributes
    ----------
    type : str
        Always ``"Cavern_T"``.
    T_values : list of float
        Time-sampled prescribed temperatures.
    time_values : list of float
        Sample times corresponding to ``T_values``.
    T : float
        Current cavern temperature.
    T_hist : list of float
        Recorded temperature history.
    t_hist : list of float
        Recorded time history.
    """

    def __init__(
        self, cavern_name: str, T_values: list, time_values: list, h_conv: float = None
    ):
        super().__init__(cavern_name, None, h_conv)
        self.type = "Cavern_T"
        self.T_values = T_values
        self.time_values = time_values
        self.T = self.T_values[0]

        # Initialize histories
        self.T_hist = [self.T]
        self.t_hist = [self.time_values[0]]

    def update_cavern(self, t: float, dt: float) -> None:
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")

    def record_data(self, t: float) -> None:
        self.T_hist.append(self.T)
        self.t_hist.append(t)

    def create_data(self, output_dir: str):
        pass
