# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Automatic differentiation backend built on `torch.func`."""

from __future__ import annotations
from typing import Callable
import torch as to

from .autodiff import ADEvaluatorBase


class TorchADEvaluator(ADEvaluatorBase):
    """
    Exact derivatives via `torch.func.jvp` (forward-mode dual numbers).

    Needs no extra dependency and evaluates the kernels in the same namespace
    (`torch`) the finite-difference backend uses, so any mechanism whose rate
    law is written out-of-place works here without further changes.
    """

    name = "torch_ad"

    def _jvp(self, fn: Callable, x: to.Tensor, v: to.Tensor) -> tuple:
        return to.func.jvp(lambda a: fn(to, a), (x,), (v,))

    def _jvp_batch(self, fn: Callable, x: to.Tensor, directions: to.Tensor) -> tuple:
        """Batched JVP using vmap over the tangent directions."""
        def single_jvp(v):
            return to.func.jvp(lambda a: fn(to, a), (x,), (v,))

        values, tangents = to.func.vmap(single_jvp)(directions)
        return values[0], tangents
