# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause
"""
SafeInCave
=========

A FEniCSx-based 3D simulator designed to simulate the mechanical behavior of salt
caverns under different operational conditions.

This module exposes the public API for the package and
sets version information.
"""

__version__ = "3.0.3"

from .Grid import GridHandlerGMSH
from .HeatEquation import HeatDiffusion
from .MomentumEquation import LinearMomentumBase, LinearMomentum, LinearMomentumMixed
from .MaterialProps import Material
from .ConstitutiveModels import *  # noqa: F403, F405
from .ConstitutiveModels import __all__ as _CONSTITUTIVE_ALL
from .OutputHandler import SaveFields
from .SimulationLogging import (
    SimulationLogging,
    register_variable,
    get_variable,
    list_registered_variables,
)
from .Thermodynamics import CavernThermodynamics
from .ConvergenceCriteria import (
    ConvergenceCriterion,
    StrainBasedCriterion,
    ForceResidualCriterion,
    DisplacementIncrementCriterion,
    ForceDisplacementCriterion,
)
from . import CavernBC
from .Simulators import (
    Simulator_TM,
    Simulator_T,
    Simulator_M,
)
from .ScreenOutput import ScreenPrinter
from .TimeHandler import (
    TimeControllerBase,
    TimeController,
    TimeControllerParabolic,
    TimeControllerAdaptive,
)
from . import MomentumBC
from . import HeatBC
from . import PostProcessingTools
from . import Utils


__all__ = [
    "GridHandlerGMSH",
    "HeatDiffusion",
    "LinearMomentumBase",
    "LinearMomentum",
    "LinearMomentumMixed",
    "Material",
    "SaveFields",
    "SimulationLogging",
    "register_variable",
    "get_variable",
    "list_registered_variables",
    "Simulator_TM",
    "Simulator_T",
    "Simulator_M",
    "ScreenPrinter",
    "TimeControllerBase",
    "TimeController",
    "TimeControllerParabolic",
    "TimeControllerAdaptive",
    "MomentumBC",
    "HeatBC",
    "CavernBC",
    "CavernThermodynamics",
    "PostProcessingTools",
    "Utils",
    "ConvergenceCriterion",
    "StrainBasedCriterion",
    "ForceResidualCriterion",
    "DisplacementIncrementCriterion",
    "ForceDisplacementCriterion",
] + _CONSTITUTIVE_ALL

__author__ = "Hermínio T. Honório"
__email__ = "h.tasinafohonorio@tno.nl"
