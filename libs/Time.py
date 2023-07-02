import numpy as np

# class TimeHandler():

# 	def __init__(self, time_settings):
# 		self.time_settings = time_settings
# 		self.build_times()

# 	def build_times(self):
# 		if type(self.time_settings["timeList"]) == type(None):
# 			self.final_time = self.time_settings["finalTime"]
# 			self.time_step = self.time_settings["timeStep"]
# 			self.time_list = np.arange(0, self.final_time + self.time_step, self.time_step)
# 		else:
# 			self.time_list = self.time_settings["timeList"]
# 			self.final_time = self.time_list[-1]
# 		self.time_step_list = np.diff(self.time_list)
# 		self.time = self.time_list[0]
# 		self.time_step = self.time_step_list[0]
# 		self.idx = 0

# 	def advance_time(self):
# 		self.idx += 1
# 		self.time = self.time_list[self.idx]
# 		self.time_step = self.time_step_list[self.idx-1]

# 	def is_final_time_reached(self):
# 		if self.time < self.final_time:
# 			return False
# 		else:
# 			return True

class TimeHandler():
	_instance = None

	def __new__(self, *args, **kwargs):
		if not self._instance:
			self._instance = super().__new__(self)
		return self._instance

	def __init__(self, time_settings):
		self.time_settings = time_settings
		self.build_times()

	def build_times(self):
		if type(self.time_settings["timeList"]) == type(None):
			self.final_time = self.time_settings["finalTime"]
			self.time_step = self.time_settings["timeStep"]
			self.time_list = np.arange(0, self.final_time + self.time_step, self.time_step)
		else:
			self.time_list = self.time_settings["timeList"]
			self.final_time = self.time_list[-1]
		self.time_step_list = np.diff(self.time_list)
		self.time = self.time_list[0]
		self.time_step = self.time_step_list[0]
		self.idx = 0

	def advance_time(self):
		self.idx += 1
		self.time = self.time_list[self.idx]
		self.time_step = self.time_step_list[self.idx-1]

	def is_final_time_reached(self):
		if self.time < self.final_time:
			return False
		else:
			return True