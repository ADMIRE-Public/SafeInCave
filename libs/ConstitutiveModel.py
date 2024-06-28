import torch as to
from Utils import *
from Elements import *

class ConstitutiveModelHandler():
	def __init__(self, grid, input_file):
		self.n_elems = grid.mesh.num_cells()
		self.input_model = input_file["material_properties"]
		self.theta = input_file["time_settings"]["theta"]

		self.elems_ie = []
		self.elems_ve = []
		self.elems_e = []

		self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		
		self.GT_ve = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_t_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.GT_ie = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.eps_t_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.C0_inv = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.C0 = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.CT = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.stress = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		self.stress_k = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.initialize()

	def initialize(self):
		# Add elastic element(s)
		elems_e = self.__get_list_of_elements(element_class="Elastic")
		for elem_e in elems_e:
			self.add_elastic_element(elem_e)

		# Add viscoelastic element
		elems_ve = self.__get_list_of_elements(element_class="Viscoelastic")
		for elem_ve in elems_ve:
			self.add_viscoelastic_element(elem_ve)

		# Add viscoelastic element
		elems_ie = self.__get_list_of_elements(element_class="Inelastic")
		for elem_ie in elems_ie:
			self.add_inelastic_element(elem_ie)

		for elem_e in self.elems_e:
			elem_e.initialize()
			self.C0_inv += elem_e.C0_inv
			self.C0 += elem_e.C0
		
		for elem_ve in self.elems_ve:
			elem_ve.initialize()

	def add_elastic_element(self, elem_e):
		self.elems_e.append(elem_e)

	def add_viscoelastic_element(self, elem_ve):
		self.elems_ve.append(elem_ve)

	def add_inelastic_element(self, elem_ie):
		self.elems_ie.append(elem_ie)

	def compute_GT_BT_ve(self, dt):
		self.GT_ve = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		phi2 = dt*(1 - self.theta)
		for elem_ve in self.elems_ve:
			elem_ve.compute_G_B(self.stress, phi2)
			self.GT_ve += elem_ve.G
			self.BT_ve += elem_ve.B

	def compute_GT_BT_ie(self, dt):
		self.GT_ie = to.zeros((self.n_elems, 6, 6), dtype=to.float64)
		self.BT_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.elems_ie:
			elem_ie.compute_G_B(self.stress, dt)
			self.GT_ie += elem_ie.G
			self.BT_ie += elem_ie.B

	def compute_CT(self, dt):
		GT = self.GT_ie + self.GT_ve
		self.CT = to.linalg.inv(self.C0_inv + dt*(1-self.theta)*GT)

	def compute_eps_ve(self, dt):
		self.eps_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.elems_ve:
			elem_ve.compute_eps_ve(self.stress, self.stress_k, dt*(1-self.theta))
			self.eps_ve += elem_ve.eps_ve

	def compute_eps_ie(self, dt):
		self.eps_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.elems_ie:
			elem_ie.compute_eps_ie(self.stress, self.stress_k, dt*(1-self.theta))
			self.eps_ie += elem_ie.eps_ie

	def compute_eps_e(self):
		self.eps_e = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_e in self.elems_e:
			elem_e.compute_eps_e(self.stress)
			self.eps_e += elem_e.eps_e

	def compute_eps_t_ie(self, dt):
		self.eps_t_ie = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ie in self.elems_ie:
			elem_ie.compute_eps_t(dt*self.theta, dt*(1 - self.theta))
			self.eps_t_ie += elem_ie.eps_t

	def compute_eps_t_ve(self, dt):
		self.eps_t_ve = to.zeros((self.n_elems, 3, 3), dtype=to.float64)
		for elem_ve in self.elems_ve:
			elem_ve.compute_eps_t(dt*self.theta, dt*(1 - self.theta))
			self.eps_t_ve += elem_ve.eps_t

	def compute_eps_t(self, dt):
		self.compute_eps_t_ie(dt)
		self.compute_eps_t_ve(dt)
		self.eps_t = self.eps_t_ve + self.eps_t_ie


	def compute_eps_rhs(self, dt, *args):
		self.compute_eps_t(dt)
		BT = self.BT_ie + self.BT_ve
		GT = self.GT_ie + self.GT_ve
		self.eps_rhs = self.eps_t - dt*(1-self.theta)*(BT + dotdot2(GT, self.stress_k))

	def compute_stress_C0(self, eps_e):
		self.stress = dotdot2(self.C0, eps_e)

	def compute_stress(self, eps_tot, dt):
		self.compute_eps_t(dt)
		GT = self.GT_ie + self.GT_ve
		BT = self.BT_ie + self.BT_ve
		self.stress = dotdot2(self.CT, eps_tot - self.eps_t + dt*(1-self.theta)*(BT + dotdot2(GT, self.stress_k)))

	def update_stress(self):
		self.stress_k = self.stress.clone()

	def compute_eps_ie_rate(self):
		for elem_ie in self.elems_ie:
			elem_ie.compute_eps_ie_rate(self.stress, return_eps_ie=False)

	def update_eps_ie_rate_old(self):
		for elem_ie in self.elems_ie:
			elem_ie.update_eps_ie_rate_old()

	def update_eps_ie_old(self):
		for elem_ie in self.elems_ie:
			elem_ie.update_eps_ie_old()

	def compute_eps_ve_rate(self, dt):
		for elem_ve in self.elems_ve:
			elem_ve.compute_eps_ve_rate(self.stress, dt*self.theta, return_eps_ve=False)

	def update_eps_ve_rate_old(self):
		for elem_ve in self.elems_ve:
			elem_ve.update_eps_ve_rate_old()

	def update_eps_ve_old(self):
		for elem_ve in self.elems_ve:
			elem_ve.update_eps_ve_old()

	def increment_internal_variables(self, dt):
		for elem_ie in self.elems_ie:
			elem_ie.increment_internal_variables(self.stress, self.stress_k, dt)

	def update_internal_variables(self):
		for elem_ie in self.elems_ie:
			elem_ie.update_internal_variables()

	def __get_list_of_elements(self, element_class="Elastic"):
		ELEMENT_DICT = {
			"Spring": Spring,
			"KelvinVoigt": Viscoelastic,
			"DislocationCreep": DislocationCreep,
			"ViscoplasticDesai": ViscoplasticDesai
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










