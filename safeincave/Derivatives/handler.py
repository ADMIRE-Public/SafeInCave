# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Factory and global default for derivative backends."""

from __future__ import annotations
from typing import Any, Callable, Optional

from .base import DerivativeEvaluator
from .finite_difference import FiniteDifferenceEvaluator
from .torch_ad import TorchADEvaluator
from .jax_ad import JaxADEvaluator

_BACKENDS: dict = {
    "finite_difference": FiniteDifferenceEvaluator,
    "automatic_differentiation_torch": TorchADEvaluator,
    "automatic_differentiation_jax": JaxADEvaluator,
}

_DEFAULT_METHOD: str = "finite_difference"


def resolve_derivative_evaluator(
    derivative_method: "str | DerivativeEvaluator | None" = None,
) -> DerivativeEvaluator:
    """
    Resolve a derivative-method selector into a concrete strategy.

    Supported names:
    - ``"finite_difference"`` (default)
    - ``"automatic_differentiation_torch"``
    - ``"automatic_differentiation_jax"``

    Parameters
    ----------
    derivative_method : str or DerivativeEvaluator or None
        Selector name, an already-built strategy (returned unchanged), or
        `None` to use the global default.

    Returns
    -------
    DerivativeEvaluator
        Concrete backend instance.
    """
    if derivative_method is None:
        derivative_method = _DEFAULT_METHOD

    if isinstance(derivative_method, DerivativeEvaluator):
        return derivative_method
    if hasattr(derivative_method, "d_tensor_d_stress") and hasattr(
        derivative_method, "d_d_isv"
    ):
        return derivative_method

    key = str(derivative_method).strip().lower()

    if key in _BACKENDS:
        return _BACKENDS[key]()

    raise ValueError(
        f"Unknown derivative_method {derivative_method!r}. Supported values: "
        "'finite_difference', 'automatic_differentiation_torch', 'automatic_differentiation_jax'."
    )


def set_default_derivative_method(
    derivative_method: "str | DerivativeEvaluator",
) -> None:
    """
    Set the derivative backend used by mechanisms built without an explicit one.

    Only affects mechanisms constructed *after* this call.

    Examples
    --------
    >>> import safeincave as sc
    >>> sc.set_default_derivative_method("automatic_differentiation_torch")
    """
    resolve_derivative_evaluator(derivative_method)  # validate eagerly

    global _DEFAULT_METHOD
    if isinstance(derivative_method, DerivativeEvaluator):
        _DEFAULT_METHOD = derivative_method
    else:
        _DEFAULT_METHOD = str(derivative_method).strip().lower()


def get_default_derivative_method() -> Any:
    """Return the current global default derivative method."""
    return _DEFAULT_METHOD


def evaluate_derivative(
    fn: Callable,
    x: Any,
    kind: str = "tensor_d_stress",
    derivative_method: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """
    Differentiate a namespace-generic kernel with the selected backend.

    Parameters
    ----------
    fn : callable
        ``fn(xp, x)`` where `xp` is the array namespace (`torch` or
        `jax.numpy`).  See :mod:`safeincave.Derivatives.namespace`.
    x : torch.Tensor
        Evaluation point: (N, 3, 3) for the stress kinds, (N,) for ``"isv"``.
    kind : {"tensor_d_stress", "scalar_d_stress", "isv"}
        Which derivative to take:

        - ``"tensor_d_stress"`` — Jacobian of an (N, 3, 3) output → (N, 6, 6)
        - ``"scalar_d_stress"`` — gradient of an (N,) output → (N, 3, 3)
        - ``"isv"`` — derivative w.r.t. a scalar internal variable
    derivative_method : str or DerivativeEvaluator, optional
        Backend selector; `None` uses the global default
        (``"finite_difference"`` unless changed).
    **kwargs
        Forwarded to the backend (`step`, `f0`).

    Returns
    -------
    torch.Tensor or tuple
        Shape depends on `kind`; see above.
    """
    evaluator = resolve_derivative_evaluator(derivative_method)

    if kind == "tensor_d_stress":
        return evaluator.d_tensor_d_stress(fn, x, **kwargs)
    if kind == "scalar_d_stress":
        return evaluator.d_scalar_d_stress(fn, x, **kwargs)
    if kind == "isv":
        return evaluator.d_d_isv(fn, x, **kwargs)

    raise ValueError(
        f"Unknown derivative kind {kind!r}. Supported values: "
        "'tensor_d_stress', 'scalar_d_stress', 'isv'."
    )
