# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP


class Cavern(ABC):
    """
    Abstract base class for cavern boundary-condition models.

    Parameters
    ----------
    cavern_name : str
        Name/label of the cavern boundary in the mesh tags.
    fluid : str, optional
        CoolProp fluid name. If provided, it is validated against CoolProp's
        ``HEOS`` backend.
    h_conv : float, optional
        Convective heat-transfer coefficient used by the heat boundary handler.

    Attributes
    ----------
    cavern_name : str
        Boundary label associated with the cavern.
    fluid : str or None
        CoolProp fluid identifier.
    h_conv : float or None
        Convective heat-transfer coefficient.
    type : str or None
        Cavern type identifier set by subclasses.
    """

    def __init__(self, cavern_name: str, fluid: str = None, h_conv: float = None):
        self.cavern_name = cavern_name
        self.check_fluid(fluid)
        self.fluid = fluid
        self.h_conv = h_conv
        self.type = None

    def check_fluid(self, fluid: str) -> None:
        if fluid is not None:
            try:
                CP.AbstractState("HEOS", fluid)
            except ValueError:
                raise ValueError(f"Fluid '{fluid}' not recognized by CoolProp.")

    @abstractmethod
    def update_cavern(self, t: float) -> None:
        pass

    @abstractmethod
    def record_data(self, t: float) -> None:
        pass

    def calculate_initial_condition(self) -> None:
        pass

    def calculate_heat(self, T: any = None, kappa: any = None) -> None:
        pass
