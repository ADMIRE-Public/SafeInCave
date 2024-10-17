"""
Material point model. This class is experimental and needs to be tested.
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
from Elements import *
from ConstitutiveModel import *
from Utils import *
import json
import os

# MPa = 1e6
# minute = 60
# hour = 60*minute
# day = 24*hour
# year = 365*day

# def read_json(file_name):
# 	with open(file_name, "r") as j_file:
# 		data = json.load(j_file)
# 	return data

# def save_json(data, file_name):
# 	with open(file_name, "w") as f:
# 	    json.dump(data, f, indent=4)

class MaterialPoint():
	def __init__(self, input_bc, input_model):
		try:
			self.theta = input_bc["Time"]["theta"]
		except:
			self.theta = 1.0
			print(f"Warning! Theta value not found in input_bc file. theta={self.theta} is assumed.")
		self.n_elems = 1
		self.m = ConstitutiveModelHandler(self.theta, self.n_elems)
		self.time = to.tensor(input_bc["Time"]["timeList"], dtype=to.float64)
		# self.time = to.linspace(input_bc["Time"]["timeList"][0], input_bc["Time"]["timeList"][-1], 100)
		self.n_steps = len(self.time)

		self.eps_tot = to.zeros((self.n_steps, 3, 3), dtype=to.float64)

		self.__create_stresses(input_bc)
		self.build_model(input_model)

	def build_model(self, input_model):
		for elem in get_list_of_elements(input_model, self.n_elems, "Elastic"):
			self.m.add_elastic_element(elem)

		for elem in get_list_of_elements(input_model, self.n_elems, "Viscoelastic"):
			self.m.add_viscoelastic_element(elem)

		for elem in get_list_of_elements(input_model, self.n_elems, "Inelastic"):
			self.m.add_inelastic_element(elem)

	def compute_strains(self):
		# Define stress
		self.m.stress = self.stress_time[0].clone()
		self.m.update_stress()

		# Compute inelastic strain rate
		self.m.compute_eps_ie_rate()
		self.m.update_eps_ie_rate_old()

		# Compute elastic strains at t=0
		self.m.compute_eps_e()
		for elem_e in self.m.elems_e:
			self.eps_tot[0] += elem_e.eps_e[0]

		# Loop over time
		for i in range(1, len(self.time)):
		# for i in range(1, 3):
			t = float(self.time[i])
			dt = t - float(self.time[i-1])

			# Update stress
			self.m.stress = self.stress_time[i].clone()
			self.m.update_stress()

			# Iterative loop settings
			tol = 1e-16
			error = 2*tol
			ite = 0
			maxiter = 20
			while error > tol and ite < maxiter:
				ite += 1

				# Compute inelastic strain rate
				self.m.compute_eps_ie_rate()

				# Compute GT matrix field
				self.m.compute_GT_BT_ie(dt)

				# Compute eps_t = eps_ie^o + phi1*eps_rate_ie^o + phi2*eps_rate_ie
				self.m.compute_eps_t_ie(dt)

				# Compute eps_ie = self.eps_t + phi2*GT_ie:(stress - stress_k) - BT_ie
				self.m.compute_eps_ie(dt)

				# Increment internal variables of inelastic elements
				self.m.increment_internal_variables(dt)

				# Compute error
				error = np.linalg.norm(self.m.BT_ie)

			# Compute viscoelastic strain rate
			self.m.compute_GT_BT_ve(dt)
			self.m.compute_eps_ve_rate(dt)

			# Compute eps_ie_t and eps_ve_t
			self.m.compute_eps_t_ie(dt)
			self.m.compute_eps_t_ve(dt)

			# Compute eps_ne = eps_ne^o + phi1*eps_ne_rate_old + phi2*eps_ne_rate + phi2*GT_ne:(stress - stress_k) - BT_ne
			self.m.compute_eps_ie(dt)
			self.m.compute_eps_ve(dt)

			# Compute eps_e = C0^(-1):sigma
			self.m.compute_eps_e()

			# Update viscoelastic strains and strain rates
			self.m.update_eps_ve_old()
			self.m.update_eps_ve_rate_old()

			# Update inelastic strain and strain rates
			self.m.update_eps_ie_old()
			self.m.update_eps_ie_rate_old()
			# try:
			# 	# print()
			# 	print("Fvp: ", float(self.m.elems_ie[0].Fvp[0]))
			# 	print("alpha: ", float(self.m.elems_ie[0].alpha[0]))
			# 	print(self.m.elems_ie[0].eps_ie_rate[0,[0,1,2],[0,1,2]])
			# 	print(self.m.elems_ie[1].eps_ie_rate[0,[0,1,2],[0,1,2]])
			# except:
			# 	pass

			# Update qsi_old, damage, etc
			self.m.update_internal_variables()

			# Compute total strain
			self.eps_tot[i] = self.m.eps_e[0] + self.m.eps_ve[0] + self.m.eps_ie[0]

			# print(float(self.time[-1])/hour, t/hour, ite, error)

	def __create_stresses(self, input_bc):
		self.stress_time = to.zeros((self.n_steps,self.n_elems,3,3), dtype=to.float64)
		for i in range(self.n_elems):
			self.stress_time[:,i,0,0] = to.tensor(input_bc["sigma_xx"], dtype=to.float64)
			self.stress_time[:,i,1,1] = to.tensor(input_bc["sigma_yy"], dtype=to.float64)
			self.stress_time[:,i,2,2] = to.tensor(input_bc["sigma_zz"], dtype=to.float64)
			self.stress_time[:,i,0,1] = self.stress_time[:,i,1,0] = to.tensor(input_bc["sigma_xy"], dtype=to.float64)
			self.stress_time[:,i,0,2] = self.stress_time[:,i,2,0] = to.tensor(input_bc["sigma_xz"], dtype=to.float64)
			self.stress_time[:,i,1,2] = self.stress_time[:,i,2,1] = to.tensor(input_bc["sigma_yz"], dtype=to.float64)




class MaterialPoint2():
	def __init__(self, input_bc, input_model):
		try:
			self.theta = input_bc["Time"]["theta"]
		except:
			self.theta = 0.5
			print(f"Warning! Theta value not found in input_bc file. theta={self.theta} is assumed.")
		self.n_elems = 1
		self.m = ConstitutiveModelHandler(self.theta, self.n_elems)
		self.time = to.tensor(input_bc["Time"]["timeList"], dtype=to.float64)
		# self.time = to.linspace(input_bc["Time"]["timeList"][0], input_bc["Time"]["timeList"][-1], 100)
		self.n_steps = len(self.time)

		self.eps_tot = to.zeros((self.n_steps, 3, 3), dtype=to.float64)

		self.__create_stresses(input_bc)
		self.build_model(input_model)

	def build_model(self, input_model):
		for elem in get_list_of_elements(input_model, self.n_elems, "Elastic"):
			self.m.add_elastic_element(elem)

		for elem in get_list_of_elements(input_model, self.n_elems, "Viscoelastic"):
			self.m.add_viscoelastic_element(elem)

		for elem in get_list_of_elements(input_model, self.n_elems, "Inelastic"):
			self.m.add_inelastic_element(elem)

	def compute_strains(self):
		# Define stress
		self.m.stress = self.stress_time[0].clone()
		self.m.update_stress()

		# Compute old non-elastic strain rates
		self.m.compute_eps_ie_rate()
		self.m.compute_eps_ve_rate(0)

		# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.m.update_eps_ie_rate_old()
		self.m.update_eps_ve_rate_old()

		# Compute elastic strains at t=0
		self.m.compute_eps_e()
		for elem_e in self.m.elems_e:
			self.eps_tot[0] += elem_e.eps_e[0]

		# Loop over time
		for i in range(1, len(self.time)):
		# for i in range(1, 3):
			t = float(self.time[i])
			dt = t - float(self.time[i-1])

			# Update stress
			self.m.stress = self.stress_time[i].clone()
			self.m.update_stress()

			# Iterative loop settings
			tol = 1e-16
			error = 2*tol
			ite = 0
			maxiter = 20
			while error > tol and ite < maxiter:
				ite += 1

				# Compute inelastic strain rate
				self.m.compute_eps_ie_rate()

				# Compute GT matrix field
				self.m.compute_GT_BT_ie(dt)

				# Compute eps_t = eps_ie^o + phi1*eps_rate_ie^o + phi2*eps_rate_ie
				self.m.compute_eps_t_ie(dt)

				# Compute eps_ie = self.eps_t + phi2*GT_ie:(stress - stress_k) - BT_ie
				self.m.compute_eps_ie(dt)

				# Increment internal variables of inelastic elements
				self.m.increment_internal_variables(dt)

				# Compute error
				error = np.linalg.norm(self.m.BT_ie)

			# Compute viscoelastic strain rate
			self.m.compute_GT_BT_ve(dt)
			self.m.compute_eps_ve_rate(dt)

			# Compute eps_ie_t and eps_ve_t
			self.m.compute_eps_t_ie(dt)
			self.m.compute_eps_t_ve(dt)

			# Compute eps_ne = eps_ne^o + phi1*eps_ne_rate_old + phi2*eps_ne_rate + phi2*GT_ne:(stress - stress_k) - BT_ne
			self.m.compute_eps_ie(dt)
			self.m.compute_eps_ve(dt)

			# Compute eps_e = C0^(-1):sigma
			self.m.compute_eps_e()

			# Update viscoelastic strains and strain rates
			self.m.update_eps_ve_old()
			self.m.update_eps_ve_rate_old()

			# Update inelastic strain and strain rates
			self.m.update_eps_ie_old()
			self.m.update_eps_ie_rate_old()
			# try:
			# 	# print()
			# 	# print(float(self.m.elems_ie[0].Fvp[0]), float(self.m.elems_ie[0].alpha[0]), float(self.m.elems_ie[0].alpha_q[0]))
			# 	print(float(self.m.elems_ie[0].Fvp[0]), float(self.m.elems_ie[0].alpha[0]))
			# 	# print("Fvp: ", float(self.m.elems_ie[0].Fvp[0]))
			# 	# print("alpha: ", float(self.m.elems_ie[0].alpha[0]))
			# 	# print(self.m.elems_ie[0].eps_ie_rate[0,[0,1,2],[0,1,2]])
			# 	# print(self.m.elems_ie[1].eps_ie_rate[0,[0,1,2],[0,1,2]])
			# except:
			# 	pass

			# Update qsi_old, damage, etc
			self.m.update_internal_variables()

			# Compute total strain
			self.eps_tot[i] = self.m.eps_e[0] + self.m.eps_ve[0] + self.m.eps_ie[0]

			# print(float(self.time[-1])/hour, t/hour, ite, error)

	def __create_stresses(self, input_bc):
		self.stress_time = to.zeros((self.n_steps,self.n_elems,3,3), dtype=to.float64)
		for i in range(self.n_elems):
			self.stress_time[:,i,0,0] = to.tensor(input_bc["sigma_xx"], dtype=to.float64)
			self.stress_time[:,i,1,1] = to.tensor(input_bc["sigma_yy"], dtype=to.float64)
			self.stress_time[:,i,2,2] = to.tensor(input_bc["sigma_zz"], dtype=to.float64)
			self.stress_time[:,i,0,1] = self.stress_time[:,i,1,0] = to.tensor(input_bc["sigma_xy"], dtype=to.float64)
			self.stress_time[:,i,0,2] = self.stress_time[:,i,2,0] = to.tensor(input_bc["sigma_xz"], dtype=to.float64)
			self.stress_time[:,i,1,2] = self.stress_time[:,i,2,1] = to.tensor(input_bc["sigma_yz"], dtype=to.float64)





