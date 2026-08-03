# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Backtracking line search for the incremental Newton solver path."""

from __future__ import annotations
from typing import Any


class BacktrackingLineSearch:
    """
    Backtracking line search on the Newton correction (Armijo-type on the
    residual norm).

    For step lengths ``α ∈ alphas`` (descending, starting at 1), evaluate
    ``u_trial = u_base − α·δu`` and accept the first α with

        ‖R(u_trial)‖ ≤ (1 − c·α)·‖R(u_base)‖.

    Backtracks are SOLVE-FREE: each trial costs one constitutive update and
    one residual assembly (the constitutive interface stages increments
    statelessly, so trial evaluations are safe). If every α fails, the
    smallest α is kept — the iteration continues and the maxiter/cutback
    machinery decides. This damps the active-plastic-set oscillation that
    destabilized acceleration schemes on the staggered path.

    The driver must skip the search on the FIRST correction of a step, where
    the Dirichlet-gap increment has to be applied in full.
    """

    def __init__(
        self,
        alphas: tuple = (1.0, 0.5, 0.25, 0.125, 0.0625),
        c: float = 1e-4,
    ):
        self.alphas = tuple(alphas)
        self.c = float(c)
        self.last_alpha = 1.0

    def search(self, eq: Any, dt: float, r_prev: float) -> tuple:
        """
        Try the correction stored in ``eq.du_sol`` from the current ``eq.u``.

        On return, ``eq.u`` holds the accepted trial state and the
        constitutive/residual fields are already evaluated there (the
        driver can reuse them as the next iteration's evaluation).

        Returns ``(alpha, r_norm, ref_norm)``.
        """
        u_base = eq.u.x.array.copy()
        r_norm = ref_norm = None
        alpha = self.alphas[0]
        for alpha in self.alphas:
            eq.u.x.array[:] = u_base - alpha * eq.du_sol.x.array
            eq.u.x.scatter_forward()
            eq.constitutive_update(dt)
            r_norm, ref_norm = eq.assemble_residual()
            if r_norm <= (1.0 - self.c * alpha) * r_prev:
                break
        self.last_alpha = alpha
        return alpha, r_norm, ref_norm
