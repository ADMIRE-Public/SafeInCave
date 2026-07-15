# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from .adaptive import TimeControllerAdaptive


class TimeControllerNewton(TimeControllerAdaptive):
    """
    Iteration-count-driven time controller for the Newton solver path
    (Abaqus-style automatic incrementation).

    Rules:

    - Non-convergence within ``maxiter``: the simulator cuts the step back
      by ``cutback_factor`` (default ×0.25) and retries.
    - Growth ×``growth_factor`` (default 1.5) after ``growth_patience``
      (default 2) consecutive converged steps that each needed at most
      ``easy_iterations`` (default 5) Newton iterations.
    - Shrink ×``shrink_factor`` when a converged step needed more than
      ``hard_iterations`` (default: ``0.8·maxiter`` when maxiter is set,
      else 10).
    - Always clamped to ``[dt_min, dt_max]``.

    There is NO accuracy-driven dt cap on this path: the backward-Euler
    return mapping is exact per increment for perfect plasticity, so the
    step size is limited only by the Newton radius (which the iteration
    count measures directly).

    Inherits the TimeControllerAdaptive constructor; ``get_next_dt`` is
    keyword-compatible with the simulator's generic call and additionally
    accepts ``n_iterations`` (advertised via ``accepts_iterations``).
    """

    accepts_iterations = True

    def __init__(
        self,
        *args,
        cutback_factor: float = 0.25,
        easy_iterations: int = 5,
        hard_iterations: int | None = None,
        growth_patience: int = 2,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.cutback_factor = float(cutback_factor)
        self.easy_iterations = int(easy_iterations)
        self.hard_iterations = hard_iterations
        self.growth_patience = int(growth_patience)
        self._easy_streak = 0

    def _hard_limit(self) -> int:
        if self.hard_iterations is not None:
            return int(self.hard_iterations)
        if self.maxiter:
            return max(int(0.8 * self.maxiter), self.easy_iterations + 1)
        return 10

    def get_next_dt(
        self,
        convergence_ratio: float | None = None,
        n_bisections: int = 0,
        converged: bool = True,
        n_iterations: int | None = None,
    ) -> float:
        """
        Suggest the next dt from the Newton iteration count. Falls back to
        the ratio-based rule of TimeControllerAdaptive when ``n_iterations``
        is not provided.
        """
        if n_iterations is None:
            return super().get_next_dt(convergence_ratio, n_bisections, converged)

        self.last_converged = converged
        self.last_n_bisections = max(0, int(n_bisections))
        ite = int(n_iterations)

        if not converged:
            self._easy_streak = 0
            dt_candidate = self.dt * self.cutback_factor
            reason = "failed_step_cutback"
        elif ite > self._hard_limit():
            self._easy_streak = 0
            dt_candidate = self.dt * self.shrink_factor
            reason = "hard_shrink"
        elif ite <= self.easy_iterations:
            self._easy_streak += 1
            if self._easy_streak >= self.growth_patience:
                dt_candidate = self.dt * self.growth_factor
                reason = "easy_grow"
            else:
                dt_candidate = self.dt
                reason = "easy_wait"
        else:
            self._easy_streak = 0
            dt_candidate = self.dt
            reason = "moderate_keep"

        self.last_dt_reason = reason
        return self._clamp_dt(dt_candidate)
