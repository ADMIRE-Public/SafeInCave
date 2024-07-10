import torch as to
from Utils import *
import Elements as element

class ConstitutiveModel():
	"""
	This is a description example.

	Attributes:
		n_elems (int): Number of grid elements.
		input_constitutive_model (dict): Dictionary containing the settings for the constitutive model.
	"""
	def __init__(self, n_elems, input_constitutive_model):
		self.n_elems = n_elems
		self.input_model = input_constitutive_model

		self.elems_ie = []
		self.elems_ve = []
		self.elems_e = []

		self.C0_inv = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.C0 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)

		# Add elastic elements
		elems_e = self.__get_list_of_elements(element_class="Elastic")
		for elem_e in elems_e:
			self.elems_e.append(elem_e)

		# Add viscoelastic elements
		elems_ve = self.__get_list_of_elements(element_class="Viscoelastic")
		for elem_ve in elems_ve:
			self.elems_ve.append(elem_ve)

		# Add inelastic elements
		elems_ie = self.__get_list_of_elements(element_class="Inelastic")
		for elem_ie in elems_ie:
			self.elems_ie.append(elem_ie)

		for elem_e in self.elems_e:
			elem_e.initialize()
			self.C0_inv += elem_e.C0_inv
			self.C0 += elem_e.C0

	def __get_list_of_elements(self, element_class="Elastic"):
		ELEMENT_DICT = {
			"Spring": element.Spring,
			"KelvinVoigt": element.Viscoelastic,
			"DislocationCreep": element.DislocationCreep,
			"ViscoplasticDesai": element.ViscoplasticDesai
		}
		list_of_elements = []
		props = self.input_model[element_class]
		for elem_name in props.keys():
			if props[elem_name]["active"] == True:
				element_parameters = props[elem_name]["parameters"]
				for param in element_parameters:
					element_parameters[param] = to.tensor(element_parameters[param])
				elem = ELEMENT_DICT[props[elem_name]["type"]](element_parameters)
				list_of_elements.append(elem)
		if element_class == "Elastic" and len(list_of_elements) == 0:
			raise Exception("Model must have at least 1 elastic element (Spring). None was given.")
		return list_of_elements










