import abc
import numpy as np
from fenics import as_backend_type

class Controller(metaclass=abc.ABCMeta):
	@abc.abstractmethod
	def __init__(self, name):
		self.name = name
		self.variable = None

	@abc.abstractmethod
	def execute(self):
		pass

	@abc.abstractmethod
	def check(self):
		pass

	@abc.abstractmethod
	def reset(self):
		pass

class IterationController(Controller):
	def __init__(self, name, max_ite=30):
		super().__init__(name)
		self.max_ite = max_ite
		self.variable = 0

	def execute(self):
		self.variable += 1

	def reset(self):
		self.variable = 0

	def check(self):
		return self.variable < self.max_ite

class ErrorController(Controller):
	def __init__(self, name, model, tol=1e-9):
		super().__init__(name)
		self.model = model
		self.tol = tol
		self.variable = 20*tol
		self.ite = 0

	def reset(self):
		self.ite = 0
		self.variable = 20*self.tol

	def execute(self):
		self.ite += 1
		self.variable = np.linalg.norm(self.model.u_k.vector() - self.model.u.vector()) / np.linalg.norm(self.model.u.vector())

	def check(self):
		if self.ite >= 2: # Iterate at least 2 times
			return self.variable > self.tol
		else:
			return True

class TimeController(Controller):
	def __init__(self, name, time_handler):
		super().__init__(name)
		self.time_handler = time_handler
		self.variable = self.time_handler.time

	def reset(self):
		self.variable = str(self.time_handler.time) + "/" + str(self.time_handler.final_time)
		# pass

	def execute(self):
		pass

	def check(self):
		return None