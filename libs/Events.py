import abc
import os

class Event(metaclass=abc.ABCMeta):
	@abc.abstractmethod
	def initialize(self):
		pass

	@abc.abstractmethod
	def execute(self):
		pass

	@abc.abstractmethod
	def finalize(self):
		pass

class ScreenOutput(Event):
	def __init__(self):
		self.controllers = []
		self.widths = []
		self.aligns = []
		self.internal_counter = 1

	def add_controller(self, controller, width=None, align=None):
		self.controllers.append(controller)

		if width == None:
			self.widths.append(len(controller.name))
		else:
			self.widths.append(width)

		if align == None or align == "center":
			self.aligns.append("^")
		elif align == "left":
			self.aligns.append("<")
		elif align == "right":
			self.aligns.append(">")
		else:
			print("Allowed values for align are center, left or right.")

	def initialize(self):
		self.__build_header()
		self.__build_bar()
		print(self.bar)
		print(self.header)
		print(self.bar)

	def execute(self):
		self.__build_row()
		print(self.row)

	def finalize(self):
		print(self.bar)

	def __build_header(self):
		n = len(self.controllers)
		if n == 0:
			self.header = "| No controllers added |"
		else:
			self.header = "|"
		for i in range(n):
			column_width = " {" + ":" + str(self.aligns[i]) + str(self.widths[i]) + "} |"
			self.header += column_width.format(self.controllers[i].name)

	def __build_bar(self):
		n = len(self.controllers)
		if n > 0:
			self.bar = "+"
			for i in range(n):
				self.bar += (self.widths[i]+2)*"-" + "+"
		else:
			self.bar = "+" + (len(self.header)-2)*"-" + "+"

	def __build_row(self):
		n = len(self.controllers)
		if n > 0:
			self.row = "|"
			for i in range(n):
				column_width = " {" + ":" + str(self.aligns[i]) + str(self.widths[i]) + "} |"
				value = str(self.controllers[i].variable)
				# value = float(value[:self.widths[i]])
				try:
					value = float(value[:self.widths[i]])
				except:
					value = value[:self.widths[i]]
				self.row += column_width.format(value)
		else:
			column_width = "|{" + ":" + "^" + str((len(self.header)-2)) + "}|"
			self.row = column_width.format(self.internal_counter)
			self.internal_counter += 1


class AverageSaver(Event):
	def __init__(self, dx, field_name, field, time_handler, output_folder):
		from ResultsHandler import TensorSaver
		self.saver = TensorSaver(field_name, dx)
		self.field = field
		self.time_handler = time_handler
		self.output_folder = output_folder

	def initialize(self):
		pass

	def execute(self):
		self.saver.record_average(self.field, self.time_handler.time)

	def finalize(self):
		self.saver.save(os.path.join(self.output_folder, "avg"))


class AverageScalerSaver(Event):
	def __init__(self, dx, field_name, field, time_handler, output_folder):
		from ResultsHandler import ScalarSaver
		self.saver = ScalarSaver(field_name, dx)
		self.field = field
		self.time_handler = time_handler
		self.output_folder = output_folder

	def initialize(self):
		pass

	def execute(self):
		self.saver.record_average(self.field, self.time_handler.time)

	def finalize(self):
		self.saver.save(os.path.join(self.output_folder, "avg"))



class VtkSaver(Event):
	def __init__(self, field_name, field, time_handler, output_folder):
		from fenics import File
		self.field = field
		self.time_handler = time_handler
		self.output_folder = os.path.join(output_folder, "vtk")
		self.vtk = File(os.path.join(self.output_folder, f"{field_name}.pvd"))

	def initialize(self):
		pass

	def execute(self):
		self.vtk << (self.field, self.time_handler.time)

	def finalize(self):
		pass


class TimeLevelCounter(Event):
	def __init__(self, time_handler, name="Time Level"):
		self.name = name
		self.counter = 0
		self.total_steps = len(time_handler.time_list)
		self.__compose_variable()

	def initialize(self):
		pass

	def execute(self):
		self.counter += 1
		self.__compose_variable()

	def finalize(self):
		pass

	def __compose_variable(self):
		self.variable = str(self.counter) + "/" + str(self.total_steps)


class TimeCounter(Event):
	def __init__(self, time_handler, name="Time", unit="hours"):
		self.name = f"{name} ({unit})"
		if unit == "seconds":
			self.unit = 1
		elif unit == "minutes":
			self.unit = 60
		elif unit == "hours":
			self.unit = 60*60
		elif unit == "days":
			self.unit = 24*60*60
		else:
			self.unit = 1
			print(f"Unit {unit} is not valid. Choose between seconds, minutes, hours or days.")
		self.time_handler = time_handler
		self.__compose_variable()

	def initialize(self):
		pass

	def execute(self):
		self.__compose_variable()

	def finalize(self):
		pass

	def __compose_variable(self):
		self.variable = str(round(self.time_handler.time/self.unit, 2))
		self.variable += "/"
		self.variable += str(round(self.time_handler.final_time/self.unit, 2))
