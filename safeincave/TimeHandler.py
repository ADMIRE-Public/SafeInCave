# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from abc import ABC, abstractmethod
from typing import Callable
from .Utils import minute, hour, day, year
import numpy as np

# Type alias
Fn = Callable[[float], float]


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

class TimeControllerAdaptive(TimeControllerBase):
    """
    Adaptive-step time controller.

    Advances the current time by a step `dt` expressed in the chosen
    unit. dt is increased if the number of iterations are small enough or decreased if too much

    Parameters
    ----------
    initial_dt : float
    Initial time-step size expressed in the units given by `time_unit`.
    max_dt : float
    Maximum time-step size expressed in the units given by `time_unit`.
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
        initial_dt: float,
        max_dt: float,
        initial_time: float,
        final_time: float,
        time_unit: str = "second",
        iterations_min: int = 5,
        iterations_max: int = 10,
        inflation: float = 2.0,
    ):
        super().__init__(initial_time, final_time, time_unit)
        self.dt = initial_dt * self.time_conversion
        self.max_dt = max_dt * self.time_conversion
        self.flag_functionOfIteration = True
        self.iterations_min = iterations_min
        self.iterations_max = iterations_max
        self.inflation = inflation

    def advance_time(self, numberIterations: int = 0) -> None:
        """
        Increment the current time by a `dt` that changes as a function of the number of iterations.

        Returns
        -------
        None
        """

        # this numberIterations == 0 is here to make timeAdaptive controller work whenever it is called as a regular non-adaptive timecontroller
        if self.step_counter == 0 or numberIterations == 0:
            pass
        elif numberIterations<=self.iterations_min:
            self.dt=self.dt*self.inflation
            if(self.dt>self.max_dt):
                self.dt = self.max_dt
        elif numberIterations>=self.iterations_max:
            self.dt=self.dt/self.inflation

        self.step_counter += 1
        self.t += self.dt


class TimeControllerParabolic(TimeControllerBase):
    """
    Nonuniform (parabolic) time controller.

    Builds a monotonically increasing list of `n_time_steps` time instants
    between `initial_time` and `final_time` by mapping an equally spaced grid
    through a parabolic function and rescaling it back to the original range.
    The step size `dt` changes over time.

    Parameters
    ----------
    n_time_steps : int
        Number of time instants (length of the time list).
    initial_time : float
        Start time expressed in the units given by `time_unit`.
    final_time : float
        Final time expressed in the units given by `time_unit`.
    time_unit : {"second", "minute", "hour", "day", "year"}, default="second"
        Unit used to interpret `initial_time` and `final_time`.

    Attributes
    ----------
    n_time_steps : int
        Number of generated time instants.
    time_list : numpy.ndarray
        Monotone array of times (in **seconds**) of length `n_time_steps`.
    dt : float
        Current time step size (in **seconds**); initialized as
        ``time_list[1] - time_list[0]``.
    step_counter : int
        Index of the most recently advanced step (starts at 0).
    """

    def __init__(
        self,
        n_time_steps: int,
        initial_time: float,
        final_time: float,
        time_unit: str = "second",
    ):
        super().__init__(initial_time, final_time, time_unit)
        self.n_time_steps = n_time_steps
        self.time_list = self.calculate_varying_times(self.fun_parabolic)
        self.dt = self.time_list[1] - self.time_list[0]
        self.step_counter = 0

    def fun_parabolic(self, t_array: np.ndarray) -> np.ndarray:
        """
        Parabolic mapping used to skew an equally spaced time grid.

        Parameters
        ----------
        t_array : numpy.ndarray
            Input array (typically an equally spaced grid).

        Returns
        -------
        numpy.ndarray
            Elementwise square of `t_array` (i.e., ``t_array**2``).
        """
        return t_array**2

    def calculate_varying_times(self, fun: Fn) -> np.ndarray:
        """
        Generate a nonuniform time grid via a monotone mapping and rescaling.

        Steps
        -----
        1. Build an equally spaced grid on ``[t_initial, t_final]``.
        2. Apply ``fun`` to this grid to skew the spacing.
        3. Linearly rescale the mapped values back to the original range so that
           the first and last times remain ``t_initial`` and ``t_final``.

        Parameters
        ----------
        fun : Callable[[numpy.ndarray], numpy.ndarray]
            Monotone mapping to skew the spacing (e.g., :meth:`fun_parabolic`).

        Returns
        -------
        numpy.ndarray
            Monotone array of times (in seconds) of length ``n_time_steps``.
        """
        t_eq = np.linspace(self.t_initial, self.t_final, self.n_time_steps)
        y = fun(t_eq)
        f_min = np.min(t_eq)
        f_max = np.max(y)
        k = (t_eq.max() - t_eq.min()) / (f_max - f_min)
        y = k * (y - f_min) + t_eq.min()
        return y

    def advance_time(self) -> None:
        """
        Advance to the next time in :attr:`time_list`.

        Increments :attr:`step_counter`, sets :attr:`t` to the corresponding
        entry of :attr:`time_list`, and updates :attr:`dt` as the difference
        between the current and previous times.

        Returns
        -------
        None

        Notes
        -----
        No bounds checking is performed. Ensure
        ``step_counter < n_time_steps - 1`` before advancing.
        """
        self.step_counter += 1
        self.t = self.time_list[self.step_counter]
        self.dt = (
            self.time_list[self.step_counter] - self.time_list[self.step_counter - 1]
        )
