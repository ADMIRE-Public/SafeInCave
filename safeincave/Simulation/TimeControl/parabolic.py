# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

from typing import Callable
import numpy as np
from .base import TimeControllerBase

# Type alias
Fn = Callable[[float], float]


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
        self.dt_used = self.dt
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
        self.dt_used = self.dt
