from fenics import *
import numpy as np
import sympy as sy
from Utils import *
import abc

class BaseElement(metaclass=abc.ABCMeta):
	@abc.abstractmethod
	def __init__(self, fem_handler):
		self.du = fem_handler.du
		self.v = fem_handler.v
		self.dx = fem_handler.dx()
		self.TS = fem_handler.TS
		self.P0 = fem_handler.P0
		self.initialize_tensors()

	def initialize_tensors(self):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.stress = local_projection(zero_tensor, self.TS)
		self.eps_tot = local_projection(zero_tensor, self.TS)

	@abc.abstractmethod
	def build_A(self, **kwargs):
		pass

	@abc.abstractmethod
	def build_b(self, **kwargs):
		pass

	# @abc.abstractmethod
	def compute_stress(self, **kwargs):
		pass

	def compute_total_strain(self, u):
		self.eps_tot.assign(local_projection(epsilon(u), self.TS))


class ElasticElement(BaseElement):
	def __init__(self, fem_handler, input_model, element_name="Spring"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.A = 0
		self.b = 0

		self.__initialize_tensors()
		self.__load_props(input_model)
		self.__initialize_constitutive_matrices()

	def build_A(self):
		a_form = inner(sigma(self.C0, epsilon(self.du)), epsilon(self.v))*self.dx
		self.A = assemble(a_form)

	def build_b(self):
		pass

	def compute_stress(self):
		stress_form = voigt2stress(dot(self.C0, strain2voigt(self.eps_e)))
		self.stress.assign(local_projection(stress_form, self.TS))

	def compute_elastic_strain(self, eps_ie=None):
		eps = self.eps_tot
		if eps_ie != None:
			eps -= eps_ie
		self.eps_e.assign(local_projection(eps, self.TS))

	def __initialize_tensors(self):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_e = local_projection(zero_tensor, self.TS)
		self.stress = local_projection(zero_tensor, self.TS)

	def __load_props(self, input_model):
		self.E = Constant(input_model["Elements"][self.element_name]["E"])
		self.nu = Constant(input_model["Elements"][self.element_name]["nu"])

	def __initialize_constitutive_matrices(self):
		self.C0_sy = constitutive_matrix_sy(self.E, self.nu)
		self.C0 = as_matrix(Constant(np.array(self.C0_sy).astype(np.float64)))


class ViscoelasticElement(BaseElement):
	def __init__(self, fem_handler, input_model, element_name="Viscoelastic"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = input_model["Time"]["theta"]
		self.A = 0
		self.b = 0

		self.__initialize_tensors()
		self.__load_props(input_model)
		self.__initialize_constitutive_matrices()

	def build_A_elastic(self):
		a_form = inner(sigma(self.C0, epsilon(self.du)), epsilon(self.v))*self.dx
		self.A_elastic = assemble(a_form)

	def build_A(self):
		a_form = inner(sigma(self.C0 - (1-self.theta)*self.C5, epsilon(self.du)), epsilon(self.v))*self.dx
		self.A = assemble(a_form)

	def build_b(self, eps_ie=None, eps_ie_old=None):
		b_form = 0
		b_form += inner(sigma(self.C4, self.eps_v_old), epsilon(self.v))*self.dx
		b_form += inner(sigma(self.theta*self.C5, self.eps_tot_old), epsilon(self.v))*self.dx
		if eps_ie != None:
			eps_ie_theta = self.theta*eps_ie_old + (1 - self.theta)*eps_ie
			b_form -= inner(sigma(self.C5, eps_ie_theta), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def compute_stress(self):
		stress_form = voigt2stress(dot(self.C0, strain2voigt(self.eps_e)))
		self.stress.assign(local_projection(stress_form, self.TS))

	def compute_elastic_strain(self, eps_ie=None):
		eps = self.eps_tot - self.eps_v
		if eps_ie != None:
			eps -= eps_ie
		self.eps_e.assign(local_projection(eps, self.TS))

	def compute_viscoelastic_strain(self, eps_ie=None, eps_ie_old=None):
		eps_tot_theta = self.theta*self.eps_tot_old + (1 - self.theta)*self.eps_tot
		form_v = dot(self.C2, strain2voigt(self.eps_v_old))
		form_v += dot(self.C3, strain2voigt(eps_tot_theta))
		if eps_ie != None:
			eps_ie_theta = self.theta*eps_ie_old + (1 - self.theta)*eps_ie
			form_v -= dot(self.C3, strain2voigt(eps_ie_theta))
		self.eps_v.assign(local_projection(voigt2stress(form_v), self.TS))

	def update(self):
		self.eps_v_old.assign(self.eps_v)
		self.eps_tot_old.assign(self.eps_tot)

	def compute_constitutive_matrices(self, dt):
		C0_bar_sy = self.__multiply(dt/self.eta, self.C0_sy)
		C1_bar_sy = self.__multiply(dt/self.eta, self.C1_sy)
		I_C0_C1_sy = self.I_sy + self.__multiply(1-self.theta, C0_bar_sy+C1_bar_sy)
		I_C0_C1_inv_sy = I_C0_C1_sy.inv()

		C2_sy = I_C0_C1_inv_sy*(self.I_sy - self.__multiply(self.theta, C0_bar_sy+C1_bar_sy))
		C3_sy = I_C0_C1_inv_sy*C0_bar_sy
		C4_sy = self.C0_sy*C2_sy
		C5_sy = self.C0_sy*C3_sy

		self.C0 = as_matrix(Constant(np.array(self.C0_sy).astype(np.float64)))
		self.C1 = as_matrix(Constant(np.array(self.C1_sy).astype(np.float64)))
		self.C2 = as_matrix(Constant(np.array(C2_sy).astype(np.float64)))
		self.C3 = as_matrix(Constant(np.array(C3_sy).astype(np.float64)))
		self.C4 = as_matrix(Constant(np.array(C4_sy).astype(np.float64)))
		self.C5 = as_matrix(Constant(np.array(C5_sy).astype(np.float64)))

	def __initialize_tensors(self,):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_e = local_projection(zero_tensor, self.TS)
		self.eps_v = local_projection(zero_tensor, self.TS)
		self.eps_v_old = local_projection(zero_tensor, self.TS)
		self.eps_tot_old = local_projection(zero_tensor, self.TS)

	def __load_props(self, input_model):
		self.E0 = Constant(input_model["Elements"]["Spring"]["E"])
		self.nu0 = Constant(input_model["Elements"]["Spring"]["nu"])
		self.E1 = Constant(input_model["Elements"]["KelvinVoigt"]["E"])
		self.nu1 = Constant(input_model["Elements"]["KelvinVoigt"]["nu"])
		self.eta = Constant(input_model["Elements"]["KelvinVoigt"]["eta"])


	def __initialize_constitutive_matrices(self):
		self.I_sy = sy.Matrix(6, 6, np.identity(6).flatten())
		self.C1_sy = constitutive_matrix_sy(self.E1, self.nu1)
		self.C1 = as_matrix(Constant(np.array(self.C1_sy).astype(np.float64)))
		self.C0_sy = constitutive_matrix_sy(self.E0, self.nu0)
		self.C0 = as_matrix(Constant(np.array(self.C0_sy).astype(np.float64)))

	def __multiply(self, a, C):
		x = sy.Symbol("x")
		return (x*C).subs(x, a)


class KelvinElement(BaseElement):
	def __init__(self, fem_handler, settings, element_name="KelvinElement"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = settings["Time"]["theta"]
		self.__load_props(settings)
		self.__initialize_constitutive_matrices()
		self.__initialize_tensors()

	def build_A(self):
		pass

	def build_b(self, C):
		b_form = inner(sigma(C, self.eps_ie), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def update(self):
		self.eps_ie_old.assign(self.eps_ie)
		self.eps_ie_rate_old.assign(self.eps_ie_rate)

	def compute_viscous_strain(self, stress, dt):
		self.__compute_viscous_strain_rate(stress)
		self.eps_ie.assign(local_projection(self.eps_ie_old + dt*(self.theta*self.eps_ie_rate_old + (1 - self.theta)*self.eps_ie_rate), self.TS))

	def __compute_viscous_strain_rate(self, stress):
		sigma_v = stress - sigma(self.C, self.eps_ie)
		self.eps_ie_rate.assign(local_projection(sigma_v/self.eta, self.TS))

	def __load_props(self, settings):
		self.eta = Constant(settings[self.element_name]["eta"])
		self.E = Constant(settings[self.element_name]["E"])
		self.nu = Constant(settings[self.element_name]["nu"])

	def __initialize_constitutive_matrices(self):
		self.C_sy = constitutive_matrix_sy(self.E, self.nu)
		self.C = as_matrix(Constant(np.array(self.C_sy).astype(np.float64)))

	def __initialize_tensors(self,):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_ie = local_projection(zero_tensor, self.TS)
		self.eps_ie_old = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate_old = local_projection(zero_tensor, self.TS)



class DashpotElement(BaseElement):
	def __init__(self, fem_handler, settings, element_name="Dashpot"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = settings["Time"]["theta"]
		self.__load_props(settings)
		self.__initialize_tensors()

	def build_A(self):
		pass

	def build_b(self, C):
		b_form = inner(sigma(C, self.eps_ie), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def update(self):
		self.eps_ie_old.assign(self.eps_ie)
		self.eps_ie_rate_old.assign(self.eps_ie_rate)

	def compute_viscous_strain(self, stress, dt):
		self.__compute_viscous_strain_rate(stress)
		self.eps_ie.assign(local_projection(self.eps_ie_old + dt*(self.theta*self.eps_ie_rate_old + (1 - self.theta)*self.eps_ie_rate), self.TS))

	def __compute_viscous_strain_rate(self, stress):
		self.eps_ie_rate.assign(local_projection((1/self.eta)*stress, self.TS))

	def __load_props(self, settings):
		self.eta = Constant(settings["Elements"][self.element_name]["eta"])

	def __initialize_tensors(self,):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_ie = local_projection(zero_tensor, self.TS)
		self.eps_ie_old = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate_old = local_projection(zero_tensor, self.TS)



class DislocationCreep(BaseElement):
	def __init__(self, fem_handler, settings, element_name="DislocationCreep"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = settings["Time"]["theta"]
		self.__load_props(settings)
		self.__initialize_tensors()

	def build_A(self):
		pass

	def build_b(self, C):
		b_form = inner(sigma(C, self.eps_ie), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def update(self):
		self.eps_ie_old.assign(self.eps_ie)
		self.eps_ie_rate_old.assign(self.eps_ie_rate)

	def compute_viscous_strain(self, stress, dt):
		self.__compute_viscous_strain_rate(stress)
		self.eps_ie.assign(local_projection(self.eps_ie_old + dt*(self.theta*self.eps_ie_rate_old + (1 - self.theta)*self.eps_ie_rate), self.TS))

	def __compute_viscous_strain_rate(self, stress):
		s = stress - (1./3)*tr(stress)*Identity(3)
		von_Mises = sqrt((3/2.)*inner(s, s))
		self.eps_ie_rate.assign(local_projection(self.B*(von_Mises**(self.n-1))*s, self.TS))

	def __load_props(self, settings):
		self.A = Constant(settings["Elements"][self.element_name]["A"])
		self.n = Constant(settings["Elements"][self.element_name]["n"])
		self.T = Constant(settings["Elements"][self.element_name]["T"])
		self.R = 8.32		# Universal gas constant
		self.Q = 51600  	# Creep activation energy, [J/mol]
		self.B = float(self.A)*np.exp(-self.Q/(self.R*float(self.T)))

	def __initialize_tensors(self,):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_ie = local_projection(zero_tensor, self.TS)
		self.eps_ie_old = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate_old = local_projection(zero_tensor, self.TS)


class PressureSolutionCreep(BaseElement):
	def __init__(self, fem_handler, settings, element_name="PressureSolutionCreep"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = settings["Time"]["theta"]
		self.__load_props(settings)
		self.__initialize_tensors()

	def build_A(self):
		pass

	def build_b(self, C):
		b_form = inner(sigma(C, self.eps_ie), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def update(self):
		self.eps_ie_old.assign(self.eps_ie)
		self.eps_ie_rate_old.assign(self.eps_ie_rate)

	def compute_viscous_strain(self, stress, dt):
		self.__compute_viscous_strain_rate(stress)
		self.eps_ie.assign(local_projection(self.eps_ie_old + dt*(self.theta*self.eps_ie_rate_old + (1 - self.theta)*self.eps_ie_rate), self.TS))

	def __compute_viscous_strain_rate(self, stress):
		s = stress - (1./3)*tr(stress)*Identity(3)
		self.eps_ie_rate.assign(local_projection(self.B*s, self.TS))

	def __load_props(self, settings):
		self.A = Constant(settings["Elements"][self.element_name]["A"])
		self.n = Constant(settings["Elements"][self.element_name]["n"])
		self.d = Constant(settings["Elements"][self.element_name]["d"])
		self.B = float(self.A)/(self.d**self.n)

	def __initialize_tensors(self,):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_ie = local_projection(zero_tensor, self.TS)
		self.eps_ie_old = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate_old = local_projection(zero_tensor, self.TS)



class Damage(BaseElement):
	def __init__(self, fem_handler, settings, element_name="Damage"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = settings["Time"]["theta"]
		self.__load_props(settings)
		self.__initialize_fields()

	def build_A(self):
		pass

	def build_b(self, C):
		b_form = inner(sigma(C, self.eps_ie), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def update(self):
		self.eps_ie_old.assign(self.eps_ie)
		self.eps_ie_rate_old.assign(self.eps_ie_rate)
		self.D_old.assign(self.D)
		self.D_rate_old.assign(self.D_rate)

	def compute_viscous_strain(self, stress, dt):
		self.__compute_viscous_strain_rate(stress, dt)
		self.eps_ie.assign(local_projection(self.eps_ie_old + (dt/day)*(self.theta*self.eps_ie_rate_old + (1 - self.theta)*self.eps_ie_rate), self.TS))

	def __compute_viscous_strain_rate(self, stress, dt):
		stress_MPa = stress/MPa
		sigma_m = (1./3)*tr(stress_MPa)
		s = stress_MPa - sigma_m*Identity(3)
		von_Mises = sqrt((3/2.)*inner(s, s))
		sigma_eq = von_Mises*((1 + self.nu + 3*(1 - 2*self.nu)*(sigma_m/von_Mises)**2)*2/3)**0.5
		self.compute_damage(sigma_eq, dt)
		aux = self.D.vector()[:]
		print(aux.min(), aux.max(), np.mean(aux))
		self.eps_ie_rate.assign(local_projection(self.A*(von_Mises**(self.n-1))*s/((1 - self.D)**self.n), self.TS))

	def compute_damage(self, sigma_eq, dt):
		self.D_rate.assign(local_projection((sigma_eq/(self.B*(1 - self.D)))**self.r, self.P0))
		self.D.assign(local_projection(self.D_old + (dt/day)*(self.theta*self.D_rate_old + (1 - self.theta)*self.D_rate), self.P0))
		aux = self.D.vector()[:]
		aux[aux > 0.8] = 0.8
		self.D.vector()[:] = aux

	def __load_props(self, settings):
		self.B = Constant(settings["Elements"][self.element_name]["B"])
		self.r = Constant(settings["Elements"][self.element_name]["r"])
		self.n = Constant(settings["Elements"][self.element_name]["n"])
		self.nu = Constant(settings["Elements"][self.element_name]["nu0"])
		self.A = Constant(settings["Elements"][self.element_name]["A"])

	def __initialize_fields(self,):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_ie = local_projection(zero_tensor, self.TS)
		self.eps_ie_old = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate_old = local_projection(zero_tensor, self.TS)
		
		self.D = local_projection(Expression("0.0", degree=0), self.P0)
		self.D_old = local_projection(Expression("0.0", degree=0), self.P0)
		self.D_rate = local_projection(Expression("0.0", degree=0), self.P0)
		self.D_rate_old = local_projection(Expression("0.0", degree=0), self.P0)
		


class ViscoplasticElement(BaseElement):
	def __init__(self, fem_handler, settings, element_name="Viscoplastic"):
		super().__init__(fem_handler)
		self.element_name = element_name
		self.theta = settings["Time"]["theta"]
		self.__load_props(settings)
		self.__initialize_fields()
		self.__initialize_yield_and_potential_functions()

	def build_A(self):
		pass

	def build_b(self, C):
		b_form = inner(sigma(C, self.eps_ie), epsilon(self.v))*self.dx
		self.b = assemble(b_form)

	def update(self):
		self.eps_ie_old.assign(self.eps_ie)
		self.eps_ie_rate_old.assign(self.eps_ie_rate)
		self.update_hardening_parameters()

	def update_hardening_parameters(self):
		self.F_vp.vector()[:] = self.Fvp_array
		self.alpha.vector()[:] = self.alpha_array
		self.alpha_q.vector()[:] = self.alpha_q_array
		self.qsi.vector()[:] = self.qsi_array
		self.qsi_v.vector()[:] = self.qsi_v_array

	def compute_viscous_strain(self, stress, dt):
		self.__compute_viscous_strain_rate(stress, dt)
		# eps_ie_old_vec = self.eps_ie_old.vector()[:]
		# eps_ie_rate_old_vec = self.eps_ie_rate_old.vector()[:]
		# eps_ie_rate_vec = self.eps_ie_rate.vector()[:]
		# eps_ie_vec = eps_ie_old_vec + dt*(self.theta*eps_ie_rate_old_vec + (1-self.theta)*eps_ie_rate_vec)
		# self.eps_ie.vector()[:] = eps_ie_vec
		self.eps_ie.assign(local_projection(self.eps_ie_old + dt*(self.theta*self.eps_ie_rate_old + (1 - self.theta)*self.eps_ie_rate), self.TS))
		# self.eps_ie.assign(local_projection(self.eps_ie_old + dt*self.eps_ie_rate, self.TS))


	def __compute_viscous_strain_rate(self, stress, dt):
		self.stress_MPa.vector()[:] = -stress.vector()[:]/MPa
		self.compute_invariants(self.stress_MPa)
		self.compute_yield_surface()
		n_elems = self.alpha_q.vector()[:].size

		# Extract arrays from element fields (DG)
		strain_rates_array = np.zeros((n_elems, 3, 3))
		# Fvp_array = self.F_vp.vector()[:]
		# alpha_array = self.alpha.vector()[:]
		# alpha_q_array = self.alpha_q.vector()[:]

		# Loop over elements
		string = []
		for e in range(n_elems):

			# Get element values
			stress_elem = self.__get_tensor_at_element(self.stress_MPa, e)
			alpha_q_elem = self.__get_scalar_at_element(self.alpha_q, e)
			alpha_elem = self.__get_scalar_at_element(self.alpha, e)
			qsi_elem = self.__get_scalar_at_element(self.qsi, e)
			qsi_v_elem = self.__get_scalar_at_element(self.qsi_v, e)
			Fvp_elem = self.compute_Fvp_at_element(stress_elem, alpha_elem)

			# Initialize viscoplastic strain rate
			strain_rate = np.zeros((3,3))

			if Fvp_elem > 0.0:
				tol = 1e-16
				error = 2*tol
				maxiter = 20
				ite = 1

				Fvp_elem_last = Fvp_elem
				alpha_last = alpha_elem
				qsi_old = qsi_elem
				qsi_v_old = qsi_v_elem

				while error > tol and ite < maxiter:
					
					# Compute flow direction
					flow_direction = self.compute_dQdS_at_element(stress_elem, alpha_q_elem)

					# Compute viscoplastic strain rate
					lmbda = self.mu_1*(Fvp_elem/self.F_0)**self.N_1
					strain_rate = -lmbda*flow_direction

					# # Compute strain rate as if it was dislocation creep
					# sigma = MPa*np.array([[stress_elem[0], stress_elem[3], stress_elem[4]],
					# 				      [stress_elem[3], stress_elem[1], stress_elem[5]],
					# 				      [stress_elem[4], stress_elem[5], stress_elem[2]]
					# ])
					# s = sigma - (1./3)*(sigma[0,0] + sigma[1,1] + sigma[2,2])*np.eye(3)
					# von_Mises = sqrt((3/2.)*double_dot(s, s))
					# A = 1.6e-40
					# R = 8.32
					# Q = 191600
					# T = 298
					# n = 5
					# B = A*np.exp(-alpha_elem*Q/R/T)
					# strain_rate = -B*(von_Mises**(n-1))*s


					# Compute qsi
					increment = double_dot(strain_rate, strain_rate)**0.5*dt
					qsi_elem = qsi_old + increment

					# Compute qsi_v
					increment_v = (strain_rate[0,0] + strain_rate[1,1] + strain_rate[2,2])*dt
					qsi_v_elem = qsi_v_old + increment_v

					# Update alpha
					alpha_elem = self.a_1 / (qsi_elem**self.eta_1)

					# Update alpha_q
					alpha_q_elem = alpha_elem + self.k_v*(self.alpha_0 - alpha_elem)*(1 - qsi_v_elem/qsi_elem)

					# Compute yield function
					Fvp_elem = self.compute_Fvp_at_element(stress_elem, alpha_elem)

					# Compute error
					error = abs(alpha_elem - alpha_last)
					alpha_last = alpha_elem
					# error = abs(Fvp_elem - Fvp_elem_last)
					# Fvp_elem_last = Fvp_elem

					# Iteration control
					ite += 1

			self.Fvp_array[e] = Fvp_elem
			self.alpha_array[e] = alpha_elem
			self.alpha_q_array[e] = alpha_q_elem
			self.qsi_array[e] = qsi_elem
			self.qsi_v_array[e] = qsi_v_elem

			strain_rates_array[e] = strain_rate

		self.eps_ie_rate.vector()[:] = strain_rates_array.flatten()

		Fvp_ind_min, Fvp_min, Fvp_ind_max, Fvp_max, n_elems, Fvp_avg = self.__compute_min_max_avg(self.Fvp_array)
		alpha_ind_min, alpha_min, alpha_ind_max, alpha_max, n_elems, alpha_avg = self.__compute_min_max_avg(self.alpha_array)
		stress_e = self.__get_tensor_at_element(self.stress_MPa, Fvp_ind_min)

		# print("| " + "(%i, %.4e) | (%.4e) | (%.4e, %.4e) |"%(Fvp_ind_min, Fvp_min, self.alpha_array[Fvp_ind_min], stress_e[0], stress_e[2]))
		print("| " + "(%.4e %.4e %.4e) | (%.4e %.4e %.4e) |"%(strain_rate[0,0], strain_rate[1,1], strain_rate[2,2], stress_e[0], stress_e[1], stress_e[2]))

		# self.update_hardening_parameters()

	def __compute_min_max_avg(self, vec):
		ind_min = np.argmin(vec)
		ind_max = np.argmax(vec)
		n_elems = len(vec)
		return ind_min, vec.min(), ind_max, vec.max(), n_elems, np.average(vec)

	def __get_tensor_at_element(self, tensor_field, elem):
		ids = [9*elem+0, 9*elem+4, 9*elem+8, 9*elem+1, 9*elem+2, 9*elem+5]
		tensor_elem = tensor_field.vector()[ids]
		tensor_elem_filtered = np.where(np.abs(tensor_elem) < 1e-1, 0, tensor_elem)
		return tensor_elem_filtered

	def __get_scalar_at_element(self, scalar_field, elem):
		return scalar_field.vector()[elem]

	def compute_yield_surface(self):
		I1_star = self.I1 + self.sigma_t
		Sr = -self.J3*np.sqrt(27)/(2*self.J2**(3/2))
		F1 = (self.gamma*I1_star**2 - self.alpha*I1_star**self.n)
		F2 = (exp(self.beta_1*I1_star) - self.beta*Sr)**self.m_v
		F = self.J2 - F1*F2
		self.F_vp.assign(local_projection(F, self.P0))

	def compute_invariants(self, stress_MPa):
		i1 = stress_MPa[0,0] + stress_MPa[1,1] + stress_MPa[2,2]
		self.I1.assign(local_projection(i1, self.P0))

		i2 = (stress_MPa[0,0]*stress_MPa[1,1] + stress_MPa[0,0]*stress_MPa[2,2] + stress_MPa[2,2]*stress_MPa[1,1]
		   - stress_MPa[0,1]**2 - stress_MPa[0,2]**2 - stress_MPa[1,2]**2)
		self.I2.assign(local_projection(i2, self.P0))

		i3 = (stress_MPa[0,0]*stress_MPa[1,1]*stress_MPa[2,2] + 2*stress_MPa[0,1]*stress_MPa[0,2]*stress_MPa[1,2]
		   - stress_MPa[0,0]*stress_MPa[1,2]**2 - stress_MPa[1,1]*stress_MPa[0,2]**2 - stress_MPa[2,2]*stress_MPa[0,1]**2)
		self.I3.assign(local_projection(i3, self.P0))

		j2 = (1/3)*i1**2 - i2
		self.J2.assign(local_projection(j2, self.P0))

		j3 = (2/27)*i1**3 - (1/3)*i1*i2 + i3
		self.J3.assign(local_projection(j3, self.P0))

	# def compute_invariants(self, stress_MPa):
	# 	i1 = stress_MPa[0,0] + stress_MPa[1,1] + stress_MPa[2,2]
	# 	self.I1.assign(local_projection(i1, self.P0))

	# 	i2 = stress_MPa[0,0]*stress_MPa[1,1] + stress_MPa[0,0]*stress_MPa[2,2] + stress_MPa[2,2]*stress_MPa[1,1]
	# 	self.I2.assign(local_projection(i2, self.P0))

	# 	i3 = stress_MPa[0,0]*stress_MPa[1,1]*stress_MPa[2,2] 
	# 	self.I3.assign(local_projection(i3, self.P0))

	# 	j2 = (1/3)*i1**2# - i2
	# 	self.J2.assign(local_projection(j2, self.P0))

	# 	j3 = (2/27)*i1**3 - (1/3)*i1*i2 #+ i3
	# 	self.J3.assign(local_projection(j3, self.P0))

	def __load_props(self, settings):
		self.F_0 = settings[self.element_name]["F_0"]
		self.mu_1 = settings[self.element_name]["mu_1"]
		self.N_1 = settings[self.element_name]["N_1"]
		self.n = settings[self.element_name]["n"]
		self.a_1 = settings[self.element_name]["a_1"]
		self.eta_1 = settings[self.element_name]["eta_1"]
		self.beta_1 = settings[self.element_name]["beta_1"]
		self.beta = settings[self.element_name]["beta"]
		self.m_v = settings[self.element_name]["m_v"]
		self.gamma = settings[self.element_name]["gamma"]
		self.alpha_0 = settings[self.element_name]["alpha_0"]
		self.k_v = settings[self.element_name]["k_v"]
		self.sigma_t = settings[self.element_name]["sigma_t"]

	def __initialize_fields(self):
		zero_tensor = Expression((("0.0","0.0","0.0"), ("0.0","0.0","0.0"), ("0.0","0.0","0.0")), degree=0)
		self.eps_ie = local_projection(zero_tensor, self.TS)
		self.eps_ie_old = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate = local_projection(zero_tensor, self.TS)
		self.eps_ie_rate_old = local_projection(zero_tensor, self.TS)
		self.stress_MPa = local_projection(zero_tensor, self.TS)
		# self.dQdS = local_projection(zero_tensor, self.TS)

		# self.alpha_field = Function(self.P0)
		# self.alpha_field.assign(self.alpha)
		
		alpha_scalar = Expression("A", A=self.alpha_0, degree=0)
		self.alpha = local_projection(alpha_scalar, self.P0)
		self.alpha_q = local_projection(alpha_scalar, self.P0)

		initial_qsi = Constant((self.a_1/self.alpha_0)**(1/self.eta_1))
		self.qsi = local_projection(Expression("X", X=initial_qsi, degree=0), self.P0)
		self.qsi_v = local_projection(Expression("X", X=initial_qsi, degree=0), self.P0)

		self.I1 = local_projection(Expression("0.0", degree=0), self.P0)
		self.I2 = local_projection(Expression("0.0", degree=0), self.P0)
		self.I3 = local_projection(Expression("0.0", degree=0), self.P0)
		self.J1 = local_projection(Expression("0.0", degree=0), self.P0)
		self.J2 = local_projection(Expression("0.0", degree=0), self.P0)
		self.J3 = local_projection(Expression("0.0", degree=0), self.P0)
		self.F_vp = local_projection(Expression("0.0", degree=0), self.P0)

		# Get number of elements
		n_elems = len(self.F_vp.vector()[:])

		self.Fvp_array = np.zeros(n_elems)
		self.alpha_array = np.zeros(n_elems)
		self.alpha_q_array = np.zeros(n_elems)
		self.qsi_array = np.zeros(n_elems)
		self.qsi_v_array = np.zeros(n_elems)

	def __initialize_yield_and_potential_functions(self):
		# Define a sympy function
		s_xx = sy.Symbol("s_xx")
		s_yy = sy.Symbol("s_yy")
		s_zz = sy.Symbol("s_zz")
		s_xy = sy.Symbol("s_xy")
		s_xz = sy.Symbol("s_xz")
		s_yz = sy.Symbol("s_yz")
		a_q = sy.Symbol("a_q")

		I1 = s_xx + s_yy + s_zz
		I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
		I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
		J1 = 0
		J2 = (1/3)*I1**2 - I2
		J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
		Sr = -(J3*np.sqrt(27))/(2*J2**1.5)

		Q1 = (-a_q*I1**self.n + self.gamma*I1**2)
		Q2 = (sy.exp(self.beta_1*I1) - self.beta*Sr)**self.m_v
		Qvp = J2 - Q1*Q2

		I_star = I1 + self.sigma_t
		F1 = (-a_q*I_star**self.n + self.gamma*I_star**2)
		F2 = (sy.exp(self.beta_1*I_star) - self.beta*Sr)**self.m_v
		Fvp = J2 - F1*F2

		variables = (s_xx, s_yy, s_zz, s_xy, s_xz, s_yz, a_q)
		self.yield_function = sy.lambdify(variables, Fvp, "numpy")
		self.dQdSxx = sy.lambdify(variables, sy.diff(Qvp, s_xx), "numpy")
		self.dQdSyy = sy.lambdify(variables, sy.diff(Qvp, s_yy), "numpy")
		self.dQdSzz = sy.lambdify(variables, sy.diff(Qvp, s_zz), "numpy")
		self.dQdSxy = sy.lambdify(variables, sy.diff(Qvp, s_xy), "numpy")
		self.dQdSxz = sy.lambdify(variables, sy.diff(Qvp, s_xz), "numpy")
		self.dQdSyz = sy.lambdify(variables, sy.diff(Qvp, s_yz), "numpy")

	def compute_Fvp_at_element(self, stress_MPa, alpha):
		s_xx = stress_MPa[0]
		s_yy = stress_MPa[1]
		s_zz = stress_MPa[2]
		s_xy = stress_MPa[3]
		s_xz = stress_MPa[4]
		s_yz = stress_MPa[5]
		I1 = s_xx + s_yy + s_zz
		I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
		J2 = (1/3)*I1**2 - I2
		if J2 == 0.0:
			return 0.0
		else:
			return self.yield_function(*stress_MPa, alpha)
		# return self.yield_function(*stress_MPa, alpha)

	def compute_dQdS_at_element(self, stress_MPa, alpha):
		dQdS = np.zeros((3,3))
		dQdS[0,0] = self.dQdSxx(*stress_MPa, alpha)
		dQdS[1,1] = self.dQdSyy(*stress_MPa, alpha)
		dQdS[2,2] = self.dQdSzz(*stress_MPa, alpha)
		dQdS[0,1] = dQdS[1,0] = 0.0*self.dQdSxy(*stress_MPa, alpha)
		dQdS[0,2] = dQdS[2,0] = 0.0*self.dQdSxz(*stress_MPa, alpha)
		dQdS[1,2] = dQdS[2,1] = 0.0*self.dQdSyz(*stress_MPa, alpha)
		return dQdS


	# def compute_dQdS_at_element(self, stress_MPa, alpha):
	# 	# Finite difference derivatives
	# 	dSigma = 1e-12
	# 	dQdS = np.zeros(6)
	# 	for i in range(6):
	# 		s = stress_MPa.copy()
	# 		s[i] += dSigma
	# 		A = self.compute_Fvp_at_element(stress_MPa, alpha)
	# 		s[i] -= 2*dSigma
	# 		B = self.compute_Fvp_at_element(s, alpha)
	# 		dQdS[i] = (B - A)/dSigma/2

	# 	# print(self.__voigt2tensor(dQdS))
	# 	return self.__voigt2tensor(dQdS)

	# def __voigt2tensor(self, s):
	# 	return np.array([
	# 			[s[0], s[3], s[4]],
	# 			[s[3], s[1], s[5]],
	# 			[s[4], s[5], s[2]],
	# 		])

