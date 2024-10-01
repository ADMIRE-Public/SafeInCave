"""
It builds the constitutive model and contains all the information about it.
"""
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

import torch as to
from Utils import *
import Elements as element

class ConstitutiveModel():
	"""
    This class is responsible for constructing the constitutive model.

    Parameters
    ----------
    n_elems : int
        Number of grid elements. This is used to create arrays of C0 and C0_inv.

    input_constitutive_model : dict
        This is a dictionary that defines the constitutive model.
    """
	def __init__(self, n_elems, input_constitutive_model):
		self.input_model = input_constitutive_model

		self._elems_ie = []
		self._elems_ve = []
		self._elems_e = []
		self._n_elems = n_elems

		self._C0_inv = to.zeros((n_elems, 6, 6), dtype=to.float64)
		self._C0 = to.zeros((n_elems, 6, 6), dtype=to.float64)

		# Add elastic elements
		for elem_e in self.get_list_of_elements(element_class="Elastic"):
			self._elems_e.append(elem_e)

		# Add viscoelastic elements
		for elem_ve in self.get_list_of_elements(element_class="Viscoelastic"):
			self._elems_ve.append(elem_ve)

		# Add inelastic elements
		for elem_ie in self.get_list_of_elements(element_class="Inelastic"):
			self._elems_ie.append(elem_ie)

		for elem_e in self._elems_e:
			elem_e.initialize()
			self._C0_inv += elem_e.C0_inv
			self._C0 += elem_e.C0

	@property
	def n_elems(self):
		"""
		int : Number of grid elements.
		"""
		return self._n_elems

	@property
	def elems_e(self):
		"""
		list : This is a Python list storing elastic elements (i.e., objects of class :class:`Elements.Spring`).
		"""
		return self._elems_e

	@property
	def elems_ve(self):
		"""
		list : This is a Python list storing viscoelastic elements (i.e., objects of class :class:`Elements.Viscoelastic`).
		"""
		return self._elems_ve

	@property
	def elems_ie(self):
		"""
		list : This is a Python list storing ineelastic elements (e.g., objects of class :class:`Elements.DislocationCreep`).
		"""
		return self._elems_ie

	@property
	def C0(self):
		"""
		torch.Tensor : This is a (n_elems, 6, 6) tensor storing the stiffness matrix associated to the linear spring for each element of the grid.
		"""
		return self._C0

	@property
	def C0_inv(self):
		"""
		torch.Tensor : This is a (n_elems, 6, 6) tensor storing the inverse of C0 for each element of the grid.
		"""
		return self._C0_inv
	



	def get_list_of_elements(self, element_class="Elastic"):
		"""
		This is an internal function responsible to build a list of elements based on the input_constitutive_model dictionary.

		Parameters
		----------
		element_class : str, default: "Elastic"
			Possible values are "Elastic", "Viscoelastic" and "Inelastic".

		Returns
		-------
		list
			Python list containing the elements of a particular class (element_class).
		"""
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










