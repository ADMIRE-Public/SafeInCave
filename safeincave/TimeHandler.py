"""
Classes responsible for time management (time stepping, controlling final time, etc)
"""
# Copyright 2025 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.
from abc import ABC, abstractmethod
from typing import Callable
from Utils import minute, hour, day, year
import numpy as np

# Type alias
Fn = Callable[[float], float]

class TimeControllerBase(ABC):
	"""
    Base class for advancing and tracking simulation time.

    Converts user-provided times into a chosen unit (seconds, minutes, hours,
    days, or years) and keeps common state used by concrete controllers.

    Parameters
    ----------
    initial_time : float
        Start time expressed in the units given by `time_unit`.
    final_time : float
        Final time expressed in the units given by `time_unit`.
    time_unit : {"second", "minute", "hour", "day", "year"}, default="second"
        Unit used to interpret `initial_time` and `final_time`.

    Attributes
    ----------
    time_unit : float
        Conversion factor to seconds for the selected unit (e.g., 60 for "minute").
    t_initial : float
        Initial time in **seconds**.
    t_final : float
        Final time in **seconds**.
    t : float
        Current time in **seconds**.
    """
	def __init__(self, initial_time: float, final_time: float, time_unit: str="second"):
		self.__decide_time_unit(time_unit)
		self.t_final = final_time*self.time_unit
		self.t_initial = initial_time*self.time_unit
		self.t = initial_time*self.time_unit

	def __decide_time_unit(self, time_unit: str) -> None:
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
		if time_unit == "second":
			self.time_unit = 1
		elif time_unit == "minute":
			self.time_unit = minute
		elif time_unit == "hour":
			self.time_unit = hour
		elif time_unit == "day":
			self.time_unit = day
		elif time_unit == "year":
			self.time_unit = year
		else:
			raise Exception(f"Time unit {time_unit} not supported.")

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
    time_step : float
        Time-step size expressed in the units given by `time_unit`.
    initial_time : float
        Start time expressed in the units given by `time_unit`.
    final_time : float
        Final time expressed in the units given by `time_unit`.
    time_unit : {"second", "minute", "hour", "day", "year"}, default="second"
        Unit used to interpret `time_step`, `initial_time`, and `final_time`.

    Attributes
    ----------
    dt : float
        Fixed time-step size in **seconds**.
    """
	def __init__(self, time_step: float, initial_time: float, final_time: float, time_unit: str="second"):
		super().__init__(initial_time, final_time, time_unit)
		self.dt = time_step*self.time_unit

	def advance_time(self) -> None:
		"""
        Increment the current time by the fixed step `dt`.

        Returns
        -------
        None
        """
		self.t += self.dt



class TimeControllerParabolic():
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

    Notes
    -----
    - This class expects attributes like :attr:`t_initial`, :attr:`t_final`,
      and :attr:`time_unit` to be available (typically provided by a base
      initializer). Ensure those are set before use.
    - Advance logic uses an internal index ``self.time_step``; initialize it
      (e.g., to ``0``) before calling :meth:`advance_time`.
    """
	def __init__(self, n_time_steps: int, initial_time: float, final_time: float, time_unit: str="second"):
		super().__init__(initial_time, final_time, time_unit)
		self.n_time_steps = n_time_steps
		self.time_list = self.calculate_varying_times(self.fun_parabolic)
		self.dt = self.time_list[1] - self.time_list[0]

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

        Steps:
        1. Build an equally spaced grid on ``[t_initial, t_final]``.
        2. Apply `fun` to this grid to skew the spacing.
        3. Linearly rescale the mapped values back to the original range so
           that the first and last times remain `t_initial` and `t_final`.

        Parameters
        ----------
        fun : Callable[[float], float]
            Monotone mapping to skew the spacing (e.g., :meth:`fun_parabolic`).

        Returns
        -------
        numpy.ndarray
            Monotone array of times (in **seconds**) of length `n_time_steps`.
        """
		t_eq = np.linspace(self.t_initial, self.t_final, self.n_time_steps)
		y = fun(t_eq)
		f_min = np.min(t_eq)
		f_max = np.max(y)
		k = (t_eq.max() - t_eq.min())/(f_max - f_min)
		y = k*(y - f_min) + t_eq.min()
		return y

	def advance_time(self) -> None:
		"""
        Advance to the next time in :attr:`time_list` and update `dt`.

        Returns
        -------
        None

        Notes
        -----
        Uses an internal counter ``self.time_step`` that should be initialized
        (e.g., ``self.time_step = 0``) before the first call. After advancing,
        sets:
        - ``self.t = time_list[time_step]``
        - ``self.dt = time_list[time_step] - time_list[time_step-1]``
        """
		self.time_step += 1
		self.t = self.time_list[self.time_step]
		self.dt = self.time_list[self.time_step] - self.time_list[self.time_step-1]

