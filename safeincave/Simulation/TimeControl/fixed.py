# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from .base import TimeControllerBase


class TimeController(TimeControllerBase):
    """
    Fixed-step time controller.

    Advances the current time by a constant step `dt` expressed in the chosen
    unit.

    Parameters
    ----------
    dt : float
    Time-step size expressed in the units given by `time_unit`.
    initial_time : float
    Start time expressed in the units given by `time_unit`.
    final_time : float
    Final time expressed in the units given by `time_unit`.
    time_unit : {"second", "minute", "hour", "day", "year"}, default="second"
    Unit used to interpret `dt`, `initial_time`, and `final_time`.

    Attributes
    ----------
    dt : float
    Fixed time-step size in **seconds**.
    """

    def __init__(
        self,
        dt: float,
        initial_time: float,
        final_time: float,
        time_unit: str = "second",
    ):
        super().__init__(initial_time, final_time, time_unit)
        self.dt_used = dt * self.time_conversion
        self.dt = dt * self.time_conversion

    def advance_time(self) -> None:
        """
        Increment the current time by the fixed step `dt`.

        Returns
        -------
        None
        """
        self.step_counter += 1
        self.t += self.dt
        self.dt_used = self.dt
