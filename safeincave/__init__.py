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

# ruff: noqa: E402
# Apply local extensions BEFORE any submodule import, so extension trees can
# shadow/add safeincave modules and every import below already resolves to them.
from . import extensions

extensions.install()

from . import BC, Cavern, PostProcessing, Utils
from .Derivatives import (
    DerivativeEvaluator,
    FiniteDifferenceEvaluator,
    JaxADEvaluator,
    TorchADEvaluator,
    evaluate_derivative,
    get_default_derivative_method,
    resolve_derivative_evaluator,
    set_default_derivative_method,
)
from .Equations.Heat import HeatDiffusion
from .Equations.Momentum import LinearMomentum, LinearMomentumBase, LinearMomentumMixed
from .Materials.Constitutive import *  # noqa: F403
from .Materials.Constitutive import __all__ as _CONSTITUTIVE_ALL
from .Materials.Material import Material
from .Mesh.Grid import GridHandlerGMSH
from .Output.SaveFields import SaveFields
from .Output.Screen import ScreenPrinter
from .Output.SimLogging import (
    SimulationLogging,
    get_variable,
    list_registered_variables,
    register_variable,
)
from .Simulation.Convergence import (
    ConvergenceCriterion,
    DisplacementIncrementCriterion,
    ForceDisplacementCriterion,
    ForceResidualCriterion,
    StrainBasedCriterion,
)
from .Simulation.Simulators import (
    Simulator_M,
    Simulator_T,
    Simulator_TM,
)
from .Simulation.TimeControl import (
    TimeController,
    TimeControllerAdaptive,
    TimeControllerBase,
    TimeControllerParabolic,
)
from .Thermo.CavernThermodynamics import CavernThermodynamics

__all__ = (
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
) + tuple(_CONSTITUTIVE_ALL)

__author__ = "Hermínio T. Honório"
__email__ = "h.tasinafohonorio@tno.nl"
