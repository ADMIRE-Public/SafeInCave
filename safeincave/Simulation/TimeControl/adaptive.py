# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from .base import TimeControllerBase


class TimeControllerAdaptive(TimeControllerBase):
    """
    Adaptive-step time controller.

    Uses the relative-convergence interface via `get_next_dt(convergence_ratio, ...)`
    to decide dt growth or shrinkage based on the ratio of solver iterations
    to the maximum allowed iterations.

    Parameters
    ----------
    initial_time : float
    Start time expressed in the units given by `time_unit`.
    initial_dt : float
    Initial time-step size expressed in the units given by `time_unit`.
    final_time : float
    Final time expressed in the units given by `time_unit`.
    time_unit : {"second", "minute", "hour", "day", "year"}, default="hour"
    Unit used to interpret all time-related parameters.
    dt_min : float, default=0.001
    Minimum time-step size expressed in the units given by `time_unit`.
    dt_max : float, default=1
    Maximum time-step size expressed in the units given by `time_unit`.
    growth_factor : float, default=1.5
    Multiplicative factor for dt growth in easy-convergence regions.
    shrink_factor : float, default=0.5
    Multiplicative factor for dt shrinking in hard-convergence regions.
    easy_ratio_threshold : float, default=0.25
    Relative-iteration threshold below which convergence is considered easy.
    hard_ratio_threshold : float, default=0.50
    Relative-iteration threshold above which convergence is considered hard.
    max_bisections : int, default=5
    Maximum retry/bisection attempts for failed steps.
    maxiter : int | None, default=None
    Optional nonlinear iteration cap consumed by simulator loops.

    Attributes
    ----------
    dt : float
        Current time-step size in **seconds**.
    """

    def __init__(
        self,
        initial_time: float,
        initial_dt: float,
        final_time: float,
        time_unit: str = "hour",
        dt_min: float = 0.001,
        dt_max: float = 1,
        shrink_factor: float = 0.5,
        growth_factor: float = 1.5,
        easy_ratio_threshold: float = 0.25,
        hard_ratio_threshold: float = 0.50,
        max_bisections: int = 5,
        maxiter: int | None = None,
    ):
        super().__init__(initial_time, final_time, time_unit)

        # Time state (stored in seconds)
        self.dt_used = initial_dt * self.time_conversion
        self.dt = initial_dt * self.time_conversion
        self.dt_min = dt_min * self.time_conversion
        self.dt_max = dt_max * self.time_conversion

        self.flag_functionOfIteration = True

        # Relative-iteration controls
        self.growth_factor = growth_factor
        self.shrink_factor = shrink_factor
        self.easy_ratio_threshold = easy_ratio_threshold
        self.hard_ratio_threshold = hard_ratio_threshold
        self.max_bisections = max_bisections
        self.maxiter = None if maxiter is None else int(maxiter)

        # Diagnostics
        self.last_ratio = 0.0
        self.last_dt_reason = "initialize"
        self.last_converged = True
        self.last_n_bisections = 0

    def _clamp_dt(self, dt_value: float) -> float:
        """Clamp dt to [dt_min, dt_max]."""
        return max(self.dt_min, min(dt_value, self.dt_max))

    def get_next_dt(
        self,
        convergence_ratio: float | None,
        n_bisections: int = 0,
        converged: bool = True,
    ) -> float:
        """
        Suggest next dt from relative convergence effort.

        Parameters
        ----------
        convergence_ratio : float | None
            Relative convergence effort r = ite / maxiter.
        n_bisections : int, default=0
            Number of retries already used in the current step.
        converged : bool, default=True
            Whether the last nonlinear step converged.

        Returns
        -------
        float
            Suggested next time step (seconds), bounded by dt_min and dt_max.
        """
        self.last_converged = converged
        self.last_n_bisections = max(0, int(n_bisections))

        ratio = 0.0 if convergence_ratio is None else max(0.0, float(convergence_ratio))
        self.last_ratio = ratio

        # Decide base candidate
        if not converged:
            dt_candidate = self.dt * self.shrink_factor
            reason = "failed_step_shrink"

        elif ratio < self.easy_ratio_threshold:
            dt_candidate = self.dt * self.growth_factor
            reason = "easy_grow"

        elif ratio > self.hard_ratio_threshold:
            dt_candidate = self.dt * self.shrink_factor
            reason = "hard_shrink"

        else:
            dt_candidate = self.dt
            reason = "moderate_keep"

        # Apply bisection penalty
        if self.last_n_bisections:
            dt_candidate *= self.shrink_factor ** self.last_n_bisections
            reason += "_bisect"

        self.last_dt_reason = reason
        return self._clamp_dt(dt_candidate)

    def advance_time(self) -> None:
        """
        Increment the current time by `dt`.

        Returns
        -------
        None
        """

        remaining_time = self.t_final - self.t
        if remaining_time <= 0.0:
            return

        self.dt = min(self.dt, remaining_time)

        self.step_counter += 1
        self.t += self.dt
        self.dt_used = self.dt
