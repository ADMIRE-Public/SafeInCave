# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Pluggable derivative evaluation for SafeInCave.

Constitutive mechanisms need the sensitivity of their strain rate (and, for
models with an internal state variable, of their hardening residue) to
stress.  Historically each mechanism computed these by hand-rolled finite
differences with hard-coded step sizes.  This subpackage extracts that into a
strategy that can be swapped for exact automatic differentiation:

- ``"finite_difference"`` — the default; same numerics as before.
- ``"torch_ad"`` — exact, via `torch.func.jvp`; no extra dependency.
- ``"jax_ad"`` — exact, via `jax.jvp`; needs the optional `jax` extra.

Select per mechanism, or globally::

    import safeincave as sc

    sc.set_default_derivative_method("torch_ad")        # whole script
    creep = sc.DislocationCreep(A, Q, n, derivative_method="jax_ad")   # one element

The AD backends evaluate *namespace-generic* kernels — see
:mod:`safeincave.Derivatives.namespace`.
"""

from .base import (  # noqa: F401
    DerivativeEvaluator,
    STRESS_BASIS,
    VOIGT_ROWS_I,
    VOIGT_ROWS_J,
)
from .finite_difference import FiniteDifferenceEvaluator  # noqa: F401
from .autodiff import ADEvaluatorBase  # noqa: F401
from .torch_ad import TorchADEvaluator  # noqa: F401
from .jax_ad import JaxADEvaluator, jax_is_available  # noqa: F401
from .handler import (  # noqa: F401
    resolve_derivative_evaluator,
    set_default_derivative_method,
    get_default_derivative_method,
    evaluate_derivative,
)
from .voigt import pack_H_voigt  # noqa: F401
from .namespace import (  # noqa: F401
    build_sym3,
    split_sym3,
    deviator,
    to_namespace,
    from_namespace,
)

__all__ = [
    "DerivativeEvaluator",
    "FiniteDifferenceEvaluator",
    "ADEvaluatorBase",
    "TorchADEvaluator",
    "JaxADEvaluator",
    "jax_is_available",
    "resolve_derivative_evaluator",
    "set_default_derivative_method",
    "get_default_derivative_method",
    "evaluate_derivative",
    "pack_H_voigt",
    "build_sym3",
    "split_sym3",
    "deviator",
    "to_namespace",
    "from_namespace",
    "STRESS_BASIS",
    "VOIGT_ROWS_I",
    "VOIGT_ROWS_J",
]
