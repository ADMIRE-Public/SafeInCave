# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Convergence criteria for nonlinear iteration loops."""

from .base import ConvergenceCriterion  # noqa: F401
from .strain_based import StrainBasedCriterion  # noqa: F401
from .force_residual import ForceResidualCriterion  # noqa: F401
from .displacement_increment import DisplacementIncrementCriterion  # noqa: F401
from .force_displacement import ForceDisplacementCriterion  # noqa: F401
from .handler import (  # noqa: F401
    ConvergenceErrorHandler,
    resolve_convergence_criterion,
    initialize_convergence_state,
    compute_error_from_criterion,
)

__all__ = [
    "ConvergenceCriterion",
    "StrainBasedCriterion",
    "ForceResidualCriterion",
    "DisplacementIncrementCriterion",
    "ForceDisplacementCriterion",
    "ConvergenceErrorHandler",
    "resolve_convergence_criterion",
    "initialize_convergence_state",
    "compute_error_from_criterion",
]
