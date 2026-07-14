# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

# Backward-compatible shim; actual implementation moved to Simulation.Convergence
from .Simulation.Convergence import (  # noqa: F401
    ConvergenceCriterion,
    StrainBasedCriterion,
    ForceResidualCriterion,
    DisplacementIncrementCriterion,
    ForceDisplacementCriterion,
    ConvergenceErrorHandler,
    resolve_convergence_criterion,
    initialize_convergence_state,
    compute_error_from_criterion,
)
