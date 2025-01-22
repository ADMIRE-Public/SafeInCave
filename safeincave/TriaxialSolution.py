import numpy as np
import pandas as pd
from Utils import save_json
from ScreenOutput import ScreenPrinter
import time
import os

sec = 1.
minute = 60*sec
hour = 60*minute
day = 24*hour
kPa = 1e3
MPa = 1e6
GPa = 1e9




def double_dot(A, B):
	n, m = A.shape
	value = 0.0
	for i in range(n):
		for j in range(m):
			value += A[i,j]*B[j,i]
	return value

def trace(s):
	return s[0,0] + s[1,1] + s[2,2]

def build_stress(s_xx=0, s_yy=0, s_zz=0, s_xy=0, s_xz=0, s_yz=0):
	return np.array([s_xx, s_yy, s_zz, s_xy, s_xz, s_yz])

def voigt2tensor(s):
	return np.array([
			[s[0], s[3], s[4]],
			[s[3], s[1], s[5]],
			[s[4], s[5], s[2]],
		])

class BaseSolution():
	def __init__(self, time_list, sigmas):
		self.sigmas = sigmas
		self.time_list = time_list

	def build_stiffness_matrix(self, E, nu):
		lame = E*nu/((1+nu)*(1-2*nu))
		G = E/(2 +2*nu)
		x = 2 # Acho que aqui deveria ser dois
		self.C = np.array([
					[2*G + lame, 	lame, 		lame, 		0.,		0., 	0.],
					[lame, 			2*G + lame, lame, 		0.,		0., 	0.],
					[lame, 			lame, 		2*G + lame, 0., 	0., 	0.],
					[0., 			0., 		0., 		x*G,	0., 	0.],
					[0., 			0., 		0., 		0., 	x*G, 	0.],
					[0., 			0., 		0., 		0., 	0.,		x*G ],
				])
		self.D = np.linalg.inv(self.C)


class Elastic(BaseSolution):
	def __init__(self, input_file, element_name, time_list, sigmas):
		super().__init__(time_list, sigmas)
		self.element_name = element_name
		self.__load_properties(input_file)
		self.build_stiffness_matrix(self.E, self.nu)

	def __load_properties(self, input_file):
		self.E = input_file["constitutive_model"]["elastic"][self.element_name]["parameters"]["E"]
		self.nu = input_file["constitutive_model"]["elastic"][self.element_name]["parameters"]["nu"]

	def compute_strains(self):
		self.eps = []
		for i in range(len(self.time_list)):
			eps_value = voigt2tensor(np.dot(self.D, self.sigmas[i]))
			self.eps.append(eps_value)
		self.eps = np.array(self.eps)


class Viscoelastic(BaseSolution):
	def __init__(self, input_file, element_name, time_list, sigmas):
		super().__init__(time_list, sigmas)
		self.element_name = element_name
		self.__load_properties(input_file)
		self.build_stiffness_matrix(self.E, self.nu)
		self.build_stress_increments()
		self.sigma_0 = self.sigmas[0]

	def __load_properties(self, input_file):
		self.E = input_file["constitutive_model"]["viscoelastic"][self.element_name]["parameters"]["E"]
		self.nu = input_file["constitutive_model"]["viscoelastic"][self.element_name]["parameters"]["nu"]
		self.eta = input_file["constitutive_model"]["viscoelastic"][self.element_name]["parameters"]["eta"]

	def A(self, t):
		return (1 - np.exp(-self.E*t/self.eta))

	def build_stress_increments(self):
		self.d_sigmas = []
		for i in range(1, len(self.sigmas)):
			self.d_sigmas.append(self.sigmas[i] - self.sigmas[i-1])
		self.d_sigmas = np.array(self.d_sigmas)

	def compute_strains(self):
		self.eps = []
		shape = self.d_sigmas[0].shape
		for i in range(0, len(self.time_list)):
			values = np.zeros(shape)
			A_list = self.A(self.time_list[i] - self.time_list[:i])
			A_list = A_list.reshape((len(A_list),1))
			values = A_list*self.d_sigmas[0:i]
			soma = np.sum(values, axis=0)
			eps_value = self.A(self.time_list[i])*voigt2tensor(np.dot(self.D, self.sigma_0))
			eps_value += voigt2tensor(np.dot(self.D, soma))
			self.eps.append(eps_value)
		self.eps = np.array(self.eps)


class DislocationCreep(BaseSolution):
	def __init__(self, input_file, element_name, time_list, sigmas):
		super().__init__(time_list, sigmas)
		self.element_name = element_name
		self.__load_properties(input_file)

		self.eps_cr = np.zeros((3,3))
		self.eps_cr_old = np.zeros((3,3))
		self.eps_cr_rate = np.zeros((3,3))

	def __load_properties(self, input_file):
		A = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["A"]
		Q = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["Q"]
		T = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["T"]
		R = 8.32
		self.B = A*np.exp(-Q/R/T)
		self.n = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["n"]

	def compute_eps_cr_rate(self, sigma):
		stress = voigt2tensor(sigma)
		s = stress - (1./3)*trace(stress)*np.eye(3)
		von_Mises = np.sqrt((3/2.)*double_dot(s, s))
		self.eps_cr_rate = self.B*(von_Mises**(self.n-1))*s

	def compute_eps_cr(self, i):
		t = self.time_list[i]
		t_old = self.time_list[i-1]
		dt = t - t_old
		self.compute_eps_cr_rate(self.sigmas[i])
		self.eps_cr = self.eps_cr_old + self.eps_cr_rate*dt
		self.eps_cr_old = self.eps_cr

	def compute_strains(self):
		self.eps = [self.eps_cr]
		for i in range(1, len(self.time_list)):
			self.compute_eps_cr(i)
			self.eps.append(self.eps_cr)
		self.eps = np.array(self.eps)



class ViscoplasticDesai(BaseSolution):
	def __init__(self, input_file, element_name, time_list, sigmas):
		super().__init__(time_list, sigmas)
		self.element_name = element_name
		self.__load_properties(input_file)
		self.__initialize_variables()
		self.qsi = 0.0

	def __load_properties(self, input_file):
		self.mu_1 = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["mu_1"]
		self.N_1 = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["N_1"]
		self.n = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["n"]
		self.a_1 = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["a_1"]
		self.eta = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["eta"]
		self.a_1 = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["a_1"]
		self.beta_1 = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["beta_1"]
		self.beta = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["beta"]
		self.m = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["m"]
		self.gamma = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["gamma"]
		self.sigma_t = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["sigma_t"]
		self.alpha_0 = input_file["constitutive_model"]["inelastic"][self.element_name]["parameters"]["alpha_0"]
		self.k_v = 0.0
		self.F0 = 1

	def compute_strains(self):
		self.eps = [np.zeros((3,3))]
		for i in range(1, len(self.time_list)):
			dt = self.time_list[i] - self.time_list[i-1]
			# print(dt)
			# dt /= day
			stress_MPa = self.sigmas[i,:].copy()/MPa
			self.compute_yield_function(stress_MPa)
			# print("Fvp:", self.Fvp)
			
			ite = 1
			if self.Fvp <= 0:
				self.eps.append(self.eps[-1])
				self.alphas.append(self.alpha)
				self.alpha_qs.append(self.alpha_q)
				self.Fvp_list.append(self.Fvp)
			else:
				tol = 1e-6
				error = 2*tol
				maxiter = 60
				alpha_last = self.alpha
				while error > tol and ite < maxiter:
					strain_rate = self.__compute_strain_rate(stress_MPa)

					increment = double_dot(strain_rate, strain_rate)**0.5*dt
					self.qsi = self.qsi_old + increment


					increment_v = (strain_rate[0,0] + strain_rate[1,1] + strain_rate[2,2])*dt
					self.qsi_v = self.qsi_v_old + increment_v

					self.__update_kv(stress_MPa)
					self.__update_alpha()
					self.__update_alpha_q()

					error = abs(self.alpha - alpha_last)
					alpha_last = self.alpha
					self.compute_yield_function(stress_MPa)

					ite += 1
					if ite >= maxiter:
						print(f"Maximum number of iterations ({maxiter}) reached during time step #{i}.")

				# string = (i, self.Fvp, strain_rate[0,0], strain_rate[1,1], strain_rate[2,2], (strain_rate[0,1]+strain_rate[0,2]+strain_rate[1,2]))
				# print("%i | %.4e | %.4e | %.4e | %.4e | %.4e"%string)
				# print(i, strain_rate[0,0], strain_rate[1,1], strain_rate[2,2], (strain_rate[0,1]+strain_rate[0,2]+strain_rate[1,2]))

				self.qsi_old = self.qsi
				self.qsi_v_old = self.qsi_v
				self.eps.append(self.eps[-1] + strain_rate*dt)
				self.alphas.append(self.alpha)
				self.alpha_qs.append(self.alpha_q)
				self.Fvp_list.append(self.Fvp)
				self.qsi_list.append(self.qsi)

			# print(self.alpha, self.Fvp, ite)

		self.eps = np.array(self.eps)
		self.alphas = np.array(self.alphas)
		self.alpha_qs = np.array(self.alpha_qs)

	def __update_kv(self, stress_MPa):
		# sigma = stress_MPa[2]
		# self.k_v = -0.00085*sigma**2 + 0.015*sigma + 0.21
		# self.k_v = 0.18
		# coeffs = [2.39027657e-06, -4.54946293e-05, -6.57580943e-04,  4.99265504e-03,  1.81960713e-01, -6.45373053e-01]
		# coeffs = [-2.25330759e-07,  7.08080098e-06,  4.63967164e-05, -2.08478762e-03, -2.79699173e-02,  8.07033586e-01, -3.33527302e+00]
		# coeffs = [ 0.00859745, -0.34279313,  3.41342767]
		# coeffs = [ 0.00878372, -0.33816767,  3.28399277]
		# coeffs = [-7.66399407e-05,  3.77666296e-03, -6.25699148e-02,  4.23032481e-01, -8.83596888e-01]
		# func = np.poly1d(coeffs)
		# self.k_v = func(sigma)
		pass

	def __update_alpha(self):
		self.alpha = self.a_1 / (self.qsi**self.eta)
		# self.alpha = self.a_1 / (self.a_1/self.alpha_0 + self.qsi**self.eta)

	def __update_alpha_q(self):
		self.alpha_q = self.alpha + self.k_v*(self.alpha_0 - self.alpha)*(1 - self.qsi_v/self.qsi)

	def __compute_strain_rate(self, stress_MPa):
		n_flow = self.evaluate_flow_direction(stress_MPa, self.alpha_q)
		# print(n_flow)
		lmbda = self.mu_1*(self.Fvp/self.F0)**self.N_1
		strain_rate = lmbda*n_flow
		return strain_rate

	def __initialize_variables(self):
		self.alpha = self.alpha_0
		self.alpha_q = self.alpha_0
		self.alphas = [self.alpha]
		self.alpha_qs = [self.alpha_q]
		self.Fvp_list = [0]
		self.qsi_list = [0]
		self.qsi_old = (self.a_1/self.alpha)**(1/self.eta)
		# self.qsi_old = 0.0
		self.qsi_v_old = self.qsi_old

	def __compute_stress_invariants(self, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz):
		I1 = s_xx + s_yy + s_zz# + self.sigma_t
		I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
		I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
		return I1, I2, I3

	def __compute_deviatoric_invariants(self, I1, I2, I3):
		J1 = np.zeros(I1.size) if type(I1) == np.ndarray else 0
		J2 = (1/3)*I1**2 - I2
		J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
		return J1, J2, J3

	def __compute_Sr(self, J2, J3):
		return -(J3*np.sqrt(27))/(2*J2**1.5)

	def compute_yield_function(self, stress_MPa):
		I1, I2, I3 = self.__compute_stress_invariants(*stress_MPa)
		J1, J2, J3 = self.__compute_deviatoric_invariants(I1, I2, I3)
		if J2 == 0.0:
			self.Fvp = -100
		else:
			Sr = self.__compute_Sr(J2, J3)
			I1_star = I1 + self.sigma_t
			F1 = (-self.alpha*I1_star**self.n + self.gamma*I1_star**2)
			F2 = (np.exp(self.beta_1*I1_star) - self.beta*Sr)**self.m
			self.Fvp = J2 - F1*F2
			# print(I1_star, I1, J2, F2, F1)
			# print(self.alpha*I1_star**self.n, self.alpha, I1_star, self.n)

	def evaluate_flow_direction(self, stress_MPa, alpha_q):
		s_xx = stress_MPa[0]
		s_yy = stress_MPa[1]
		s_zz = stress_MPa[2]
		s_xy = stress_MPa[3]
		s_xz = stress_MPa[4]
		s_yz = stress_MPa[5]

		I1, I2, I3 = self.__compute_stress_invariants(*stress_MPa)
		J1, J2, J3 = self.__compute_deviatoric_invariants(I1, I2, I3)

		I1_star = I1 + self.sigma_t
		Sr = self.__compute_Sr(J2, J3)

		F1 = (-self.alpha_q*I1**self.n + self.gamma*I1**2)
		F2 = (np.exp(self.beta_1*I1) - self.beta*Sr)
		F1_star = (-self.alpha_q*I1_star**self.n + self.gamma*I1_star**2)
		F2_star = (np.exp(self.beta_1*I1_star) - self.beta*Sr)
		F = J2 - F1_star*F2_star**self.m

		I1_aux = I1#_star

		dF1_dI1 = 2*self.gamma*I1_aux - self.alpha_q*self.n*I1_aux**(self.n-1)
		dF2m_dI1 = self.beta_1*self.m*np.exp(self.beta_1*I1_aux)*F2**(self.m-1)
		dF_dI1 = -(dF1_dI1*F2**self.m + F1*dF2m_dI1)
		dF_dJ2 = 1 - (3*np.sqrt(27)*self.beta*J3*self.m*F1*F2**(self.m-1))/(4*J2**(5/2))

		dF2_dJ2 = -(3*self.beta*J3*27**0.5)/(4*J2**(5/2))
		dF_dJ2 = 1 - F1*self.m*F2**(self.m-1)*dF2_dJ2
		dF_dJ3 = -self.m*F1*self.beta*np.sqrt(27)*F2**(self.m-1)/(2*J2**1.5)

		dI1_dSxx = 1.0
		dI1_dSyy = 1.0
		dI1_dSzz = 1.0
		dI1_dSxy = 0.0
		dI1_dSxz = 0.0
		dI1_dSyz = 0.0

		dI2_dSxx = s_yy + s_zz
		dI2_dSyy = s_xx + s_zz
		dI2_dSzz = s_xx + s_yy
		dI2_dSxy = -2*s_xy
		dI2_dSxz = -2*s_xz
		dI2_dSyz = -2*s_yz

		dI3_dSxx = s_yy*s_zz - s_yz**2
		dI3_dSyy = s_xx*s_zz - s_xz**2
		dI3_dSzz = s_xx*s_yy - s_xy**2
		dI3_dSxy = 2*(s_xz*s_yz - s_zz*s_xy)
		dI3_dSxz = 2*(s_xy*s_yz - s_yy*s_xz)
		dI3_dSyz = 2*(s_xz*s_xy - s_xx*s_yz)

		dJ2_dI1 = (2/3)*I1
		dJ2_dI2 = -1.0

		dJ2_dSxx = dJ2_dI1*dI1_dSxx + dJ2_dI2*dI2_dSxx
		dJ2_dSyy = dJ2_dI1*dI1_dSyy + dJ2_dI2*dI2_dSyy
		dJ2_dSzz = dJ2_dI1*dI1_dSzz + dJ2_dI2*dI2_dSzz
		dJ2_dSxy = dJ2_dI1*dI1_dSxy + dJ2_dI2*dI2_dSxy
		dJ2_dSxz = dJ2_dI1*dI1_dSxz + dJ2_dI2*dI2_dSxz
		dJ2_dSyz = dJ2_dI1*dI1_dSyz + dJ2_dI2*dI2_dSyz

		dJ3_dI1 = (2/9)*I1**2 - (1/3)*I2
		dJ3_dI2 = -(1/3)*I1
		dJ3_dI3 = 1.0

		dJ3_dSxx = dJ3_dI1*dI1_dSxx + dJ3_dI2*dI2_dSxx + dJ3_dI3*dI3_dSxx
		dJ3_dSyy = dJ3_dI1*dI1_dSyy + dJ3_dI2*dI2_dSyy + dJ3_dI3*dI3_dSyy
		dJ3_dSzz = dJ3_dI1*dI1_dSzz + dJ3_dI2*dI2_dSzz + dJ3_dI3*dI3_dSzz
		dJ3_dSxy = dJ3_dI1*dI1_dSxy + dJ3_dI2*dI2_dSxy + dJ3_dI3*dI3_dSxy
		dJ3_dSxz = dJ3_dI1*dI1_dSxz + dJ3_dI2*dI2_dSxz + dJ3_dI3*dI3_dSxz
		dJ3_dSyz = dJ3_dI1*dI1_dSyz + dJ3_dI2*dI2_dSyz + dJ3_dI3*dI3_dSyz

		dQdS = np.zeros((3,3))
		dQdS[0,0] = dF_dI1*dI1_dSxx + dF_dJ2*dJ2_dSxx + dF_dJ3*dJ3_dSxx
		dQdS[1,1] = dF_dI1*dI1_dSyy + dF_dJ2*dJ2_dSyy + dF_dJ3*dJ3_dSyy
		dQdS[2,2] = dF_dI1*dI1_dSzz + dF_dJ2*dJ2_dSzz + dF_dJ3*dJ3_dSzz
		dQdS[0,1] = dQdS[1,0] = dF_dI1*dI1_dSxy + dF_dJ2*dJ2_dSxy + dF_dJ3*dJ3_dSxy
		dQdS[0,2] = dQdS[2,0] = dF_dI1*dI1_dSxz + dF_dJ2*dJ2_dSxz + dF_dJ3*dJ3_dSxz
		dQdS[1,2] = dQdS[2,1] = dF_dI1*dI1_dSyz + dF_dJ2*dJ2_dSyz + dF_dJ3*dJ3_dSyz

		return dQdS



model_elem_dict = {
	"Spring": Elastic,
	"DislocationCreep": DislocationCreep,
	"KelvinVoigt": Viscoelastic,
	"ViscoplasticDesai": ViscoplasticDesai
}

class MaterialPointModel:
	def __init__(self, input_file, bName_sxx, bName_syy, bName_szz):

		self.results_folder = os.path.join(input_file["output"]["path"], "material_point")
		os.makedirs(self.results_folder, exist_ok=True)

		# Screen info
		ScreenPrinter.reset_instance()
		self.screen = ScreenPrinter()
		self.screen.start_timer()
		self.screen.set_header_columns(["Running model", "Simulation time (s)"], "center")
		self.screen.set_row_formats(["%s", "%.5f"], ["center", "center"])
		self.screen.print_welcome()
		self.screen.print_comment(" ")
		self.screen.print_comment(" Material point (MP) model for triaxial test.")
		self.screen.print_comment(" ")
		self.screen.print_comment(" Results folder:")
		self.screen.print_comment(f"          {self.results_folder}")
		self.screen.print_comment(" ")

		self.load_time_list(input_file)
		self.build_sigmas(input_file, bName_sxx, bName_syy, bName_szz)

		self.screen.print_comment(" Constitutive model:")
		self.build_model(input_file, bName_sxx, bName_syy, bName_szz)

		self.screen.print_header()

	def build_model(self, input_file, bName_sxx, bName_syy, bName_szz):
		self.model_element_names = []
		self.model_elements = []
		for element_group in input_file["constitutive_model"].keys():
			for element_name in input_file["constitutive_model"][element_group].keys():
				if input_file["constitutive_model"][element_group][element_name]["active"] == True:
					element_type = input_file["constitutive_model"][element_group][element_name]["type"]
					model_object = model_elem_dict[element_type]
					self.model_elements.append(model_object(input_file, element_name, self.time_list, self.sigmas))
					self.model_element_names.append(element_name)
					self.screen.print_comment(f"          {element_name}")
		self.screen.print_comment(" ")

	def run(self):
		# Compute total strain
		self.eps_tot = 0
		self.results = {}
		for name, model_element in zip(self.model_element_names, self.model_elements):
			start = time.time()
			model_element.compute_strains()
			self.results[name] = {
				"epsilon_3": list(model_element.eps[:,0,0].copy()),
				"epsilon_1": list(model_element.eps[:,2,2].copy())
			}
			self.eps_tot += model_element.eps.copy()
			final = time.time()
			self.screen.print_row([name, final-start])
		self.results["total"] = {
			"epsilon_3": list(self.eps_tot[:,0,0]),
			"epsilon_1": list(self.eps_tot[:,2,2])
		}
		self.results["stress"] = {
			"sigma_3": list(self.sigmas[:,0]),
			"sigma_1": list(self.sigmas[:,2])
		}
		self.results["time"] = list(self.time_list)
		self.save_solution()

		self.screen.close()

		self.screen.save_log(self.results_folder)

	def save_solution(self):
		file_name = os.path.join(self.results_folder, "results.json")
		save_json(self.results, file_name)

	def load_time_list(self, input_file):
		t_final = input_file["time_settings"]["time_list"][-1]
		dt = input_file["simulation_settings"]["operation"]["dt_max"]
		self.time_list = []
		t = dt
		while t < t_final:
			self.time_list.append(t)
			t += dt
		self.time_list = np.array(self.time_list)

	def build_sigmas(self, input_file, bName_sxx, bName_syy, bName_szz):
		t_list = input_file["time_settings"]["time_list"]
		n = len(self.time_list)
		sxx_bc = input_file["boundary_conditions"][bName_sxx]["values"]
		syy_bc = input_file["boundary_conditions"][bName_syy]["values"]
		szz_bc = input_file["boundary_conditions"][bName_szz]["values"]
		sigma_xx, sigma_yy, sigma_zz = [], [], []
		for t in self.time_list:
			sigma_xx.append(np.interp(t, t_list, sxx_bc))
			sigma_yy.append(np.interp(t, t_list, syy_bc))
			sigma_zz.append(np.interp(t, t_list, szz_bc))
		sigma_xx = np.array(sigma_xx).reshape((1, n))
		sigma_yy = np.array(sigma_yy).reshape((1, n))
		sigma_zz = np.array(sigma_zz).reshape((1, n))
		sigma_xy = np.zeros_like(sigma_xx)
		sigma_yz = np.zeros_like(sigma_xx)
		sigma_xz = np.zeros_like(sigma_xx)
		self.sigmas = np.concatenate((sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz)).T
