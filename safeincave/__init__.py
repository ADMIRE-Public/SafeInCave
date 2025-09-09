"""
SafeInCave
=========

A FEniCSx-based 3D simulator designed to simulate the mechanical behavior of salt caverns under different operational conditions.

This module exposes the public API for the package and
sets version information.
"""

# Version info
__version__ = "2.0.0"

from .Grid import GridHandlerGMSH
from .HeatEquation import HeatDiffusion
from .MomentumEquation import LinearMomentumBase, LinearMomentum
from .MaterialProps import Material, NonElasticElement, Spring, Thermoelastic, Viscoelastic, DislocationCreep, PressureSolutionCreep, ViscoplasticDesai
from .OutputHandler import SaveFields
from .Simulators import Simulator_TM, Simulator_T, Simulator_M
from .ScreenOutput import ScreenPrinter
from .TimeHandler import TimeControllerBase, TimeController, TimeControllerParabolic
from . import MomentumBC
from . import HeatBC
from . import PostProcessingTools
from . import Utils


__all__ = [
    "GridHandlerGMSH",
    "HeatBC",
    "HeatDiffusion",
    "LinearMomentumBase",
    "LinearMomentum",
    "Material",
    "NonElasticElement",
    "Spring",
    "Thermoelastic",
    "Viscoelastic",
    "DislocationCreep",
    "PressureSolutionCreep",
    "ViscoplasticDesai",
    "SaveFields",
    "Simulator_TM",
    "Simulator_T",
    "Simulator_M",
    "ScreenPrinter",
    "TimeControllerBase",
    "TimeController",
    "TimeControllerParabolic",
    "MomentumBC",
    "HeatBC",
    "PostProcessingTools",
    "Utils",
]

__author__ = "Hermínio T. Honório"
__email__ = "h.tasinafohonor@tudelft.nl"

