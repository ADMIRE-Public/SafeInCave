# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
from ...Utils import minute, hour, day, year


class TimeControllerBase(ABC):
    """
    Base class for time-stepping controllers.

    Converts user-provided time units (e.g., minutes, hours) into **seconds**
    for internal bookkeeping, stores the current time, and provides a
    loop-termination predicate. Concrete subclasses implement how time
    advances each step.

    Parameters
    ----------
    initial_time : float
        Start time expressed in the unit given by ``time_unit``.
    final_time : float
        End time expressed in the unit given by ``time_unit``.
    time_unit : {"second", "minute", "hour", "day", "year"}, default="second"
        Unit of the input times. Internally, times are stored in seconds.

    Attributes
    ----------
    time_unit : str
        The time unit specified by the user.
    time_conversion : float
        Multiplicative factor to convert the chosen unit into seconds.
    t_initial : float
        Initial time in seconds.
    t_final : float
        Final time in seconds.
    t : float
        Current time in seconds.
    step_counter : int
        Number of completed steps (starts at 0).
    flag_functionOfIteration: bool
        Flag to define if the dt is a function of iteration number (adaptive)

    """

    def __init__(
        self, initial_time: float, final_time: float, time_unit: str = "second"
    ):
        self.time_unit = time_unit
        self.__decide_time_unit()
        self.t_final = final_time * self.time_conversion
        self.t_initial = initial_time * self.time_conversion
        self.t = initial_time * self.time_conversion
        self.dt = 0.0
        self.dt_used = 0.0
        self.step_counter = 0
        self.flag_functionOfIteration = False

    def __decide_time_unit(self) -> None:
        """
        Map a textual time unit to its conversion factor to seconds.

        Parameters
        ----------
        time_unit : {"second", "minute", "hour", "day", "year"}
            Name of the time unit.

        Returns
        -------
        None

        Raises
        ------
        Exception
            If `time_unit` is not one of the supported values.

        Notes
        -----
        This is an internal helper that sets :attr:`time_unit` to a scalar factor.
        """
        if self.time_unit == "second":
            self.time_conversion = 1
        elif self.time_unit == "minute":
            self.time_conversion = minute
        elif self.time_unit == "hour":
            self.time_conversion = hour
        elif self.time_unit == "day":
            self.time_conversion = day
        elif self.time_unit == "year":
            self.time_conversion = year
        else:
            raise Exception(f"Time unit {self.time_unit} not supported.")

    def keep_looping(self) -> None:
        """
        Check whether the controller should continue advancing time.

        Returns
        -------
        bool
            ``True`` while the current time `t` is strictly less than `t_final`,
            otherwise ``False``.
        """
        return self.t < self.t_final

    @abstractmethod
    def advance_time(self) -> None:
        """
        Advance the internal time by one step.

        Returns
        -------
        None

        Notes
        -----
        Subclasses must implement the time-update rule (e.g., add a fixed `dt`
        or follow a varying schedule).
        """
        pass
