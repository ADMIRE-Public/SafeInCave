# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

import os
from ..Utils import save_json
from .base import Cavern


class CavernHandler:
    """
    Container and dispatcher for cavern boundary-condition models.

    This class stores caverns by model type and forwards volume, heat,
    time-update, recording, and export operations to the registered caverns.

    Attributes
    ----------
    caverns_T : list[Cavern_T]
        Caverns with prescribed time-dependent temperature.
    caverns_PT : list[Cavern_PT]
        Caverns with prescribed time-dependent pressure and temperature.
    caverns_MFlux : list[Cavern_MassFlux]
        Caverns advanced from a prescribed mass-flux history.
    output_folder : str or None
        Folder where cavern history data are written by :meth:`save_caverns_data`.
    """

    def __init__(self):
        self.caverns_T = []
        self.caverns_PT = []
        self.caverns_MFlux = []
        self.output_folder = None

    def add_cavern(self, cavern: Cavern) -> None:
        if cavern.type == "Cavern_T":
            self.caverns_T.append(cavern)
        elif cavern.type == "Cavern_PT":
            self.caverns_PT.append(cavern)
        elif cavern.type == "Cavern_MFlux":
            self.caverns_MFlux.append(cavern)
        else:
            raise ValueError(f"Cavern type {cavern.type} not supported")

    def set_output_folder(self, output_folder: str) -> None:
        self.output_folder = output_folder

    def calculate_volumes(self, u: any = None) -> None:
        for cavern in self.caverns_MFlux + self.caverns_PT:
            cavern.calculate_volume(u)

    def calculate_total_heat(self, dt: float, T: any = None, kappa: any = None) -> None:
        for cavern in self.caverns_MFlux + self.caverns_PT:
            cavern.calculate_heat(dt, T, kappa)

    def update_caverns(
        self, t: float, dt: float = None, u: any = None, T: any = None
    ) -> None:
        for cavern in self.caverns_T:
            cavern.update_cavern(t, dt)

        for cavern in self.caverns_PT:
            cavern.update_cavern(t, dt)

        for cavern in self.caverns_MFlux:
            cavern.update_cavern(t, dt)

    def calculate_initial_conditions(self) -> None:
        for cavern in self.caverns_MFlux:
            cavern.calculate_initial_condition()

    def record_cavern_data(self, t: float) -> None:
        for cavern in self.caverns_T:
            cavern.record_data(t)

        for cavern in self.caverns_PT:
            cavern.record_data(t)

        for cavern in self.caverns_MFlux:
            cavern.record_data(t)

    def save_caverns_data(self):
        cavern_data = {}
        for cavern in self.caverns_T + self.caverns_PT + self.caverns_MFlux:
            cavern_data[cavern.cavern_name] = {}
            cavern_data[cavern.cavern_name]["Type"] = cavern.type
            cavern_data[cavern.cavern_name]["Fluid"] = cavern.fluid
            cavern_data[cavern.cavern_name]["Data"] = cavern.create_data()
        if self.output_folder is not None:
            save_json(cavern_data, os.path.join(self.output_folder, "cavern_data.json"))
