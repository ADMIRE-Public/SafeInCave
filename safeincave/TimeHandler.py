# Copyright 2024 The safeincave community.
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
from Utils import minute, hour, day, year
import numpy as np

class TimeController():
	def __init__(self, time_step, final_time, initial_time, time_unit="second"):
		self.__decide_time_unit(time_unit)
		self.dt = time_step*self.time_unit
		self.t_final = final_time*self.time_unit
		self.t = initial_time*self.time_unit

	def __decide_time_unit(self, time_unit):
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

	def advance_time(self):
		self.t += self.dt

	def keep_looping(self):
		return self.t < self.t_final


class TimeControllerParabolic():
	def __init__(self, final_time, initial_time, n_time_steps, time_unit="second"):
		self.__decide_time_unit(time_unit)
		self.n_time_steps = n_time_steps
		self.t_initial = initial_time*self.time_unit
		self.t_final = final_time*self.time_unit

		self.time_list = self.calculate_varying_times(self.fun_parabolic)
		self.time_step = 0
		self.t = self.time_list[self.time_step]
		self.dt = 0

	def __decide_time_unit(self, time_unit):
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


	def fun_parabolic(self, t_array):
		return t_array**2

	def calculate_varying_times(self, fun):
		t_eq = np.linspace(self.t_initial, self.t_final, self.n_time_steps)
		y = fun(t_eq)
		f_min = np.min(t_eq)
		f_max = np.max(y)
		k = (t_eq.max() - t_eq.min())/(f_max - f_min)
		y = k*(y - f_min) + t_eq.min()
		return y

	def advance_time(self):
		self.time_step += 1
		self.t = self.time_list[self.time_step]
		self.dt = self.time_list[self.time_step] - self.time_list[self.time_step-1]

	def keep_looping(self):
		return self.t < self.t_final


