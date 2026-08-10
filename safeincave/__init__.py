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

# Apply local extensions BEFORE any submodule import, so extension trees can
# shadow/add safeincave modules and every import below already resolves to them.
# The late imports are intentional, hence the file-level E402 exemption.
# ruff: noqa: E402
from . import extensions

extensions.install()

from .Derivatives import (
    DerivativeEvaluator,
    FiniteDifferenceEvaluator,
    TorchADEvaluator,
    JaxADEvaluator,
    resolve_derivative_evaluator,
    set_default_derivative_method,
    get_default_derivative_method,
    evaluate_derivative,
)
from .Mesh.Grid import GridHandlerGMSH, GridHandlerPythonScript
from .Equations.Heat import HeatDiffusion
from .Equations.Momentum import LinearMomentumBase, LinearMomentum, LinearMomentumMixed
from .Materials.Material import Material
from .Materials.Constitutive import *  # noqa: F403, F405
from .Materials.Constitutive import __all__ as _CONSTITUTIVE_ALL
from .Output.SaveFields import SaveFields
from .Output.DataExtract import (
    SimulationLogging,
    CompositeLogger,
    register_variable,
    get_variable,
    list_registered_variables,
    extract_point,
    extract_variable,
)
from .Thermo.CavernThermodynamics import CavernThermodynamics
from .Simulation.Convergence import (
    ConvergenceCriterion,
    StrainBasedCriterion,
    ForceResidualCriterion,
    DisplacementIncrementCriterion,
    ForceDisplacementCriterion,
)
from . import Cavern
from .Simulation.Simulators import (
    Simulator_TM,
    Simulator_T,
    Simulator_M,
    GeostaticStep,
    GeostaticReport,
    check_equilibrium,
)
from .Output.Screen import ScreenPrinter
from .Simulation.TimeControl import (
    TimeControllerBase,
    TimeController,
    TimeControllerParabolic,
    TimeControllerAdaptive,
)
from . import BC
from . import PostProcessing
from . import Utils


__all__ = [
    "DerivativeEvaluator",
    "FiniteDifferenceEvaluator",
    "TorchADEvaluator",
    "JaxADEvaluator",
    "resolve_derivative_evaluator",
    "set_default_derivative_method",
    "get_default_derivative_method",
    "evaluate_derivative",
    "GridHandlerGMSH",
    "HeatDiffusion",
    "LinearMomentumBase",
    "LinearMomentum",
    "LinearMomentumMixed",
    "Material",
    "SaveFields",
    "SimulationLogging",
    "CompositeLogger",
    "register_variable",
    "get_variable",
    "list_registered_variables",
    "extract_point",
    "extract_variable",
    "Simulator_TM",
    "Simulator_T",
    "Simulator_M",
    "GeostaticStep",
    "GeostaticReport",
    "check_equilibrium",
    "ScreenPrinter",
    "TimeControllerBase",
    "TimeController",
    "TimeControllerParabolic",
    "TimeControllerAdaptive",
    "BC",
    "Cavern",
    "CavernThermodynamics",
    "PostProcessing",
    "Utils",
    "ConvergenceCriterion",
    "StrainBasedCriterion",
    "ForceResidualCriterion",
    "DisplacementIncrementCriterion",
    "ForceDisplacementCriterion",
] + _CONSTITUTIVE_ALL

__author__ = "Hermínio T. Honório"
__email__ = "h.tasinafohonorio@tno.nl"
