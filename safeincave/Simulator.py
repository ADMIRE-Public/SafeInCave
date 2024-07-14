from ConstitutiveModel import ConstitutiveModel
from Equations import LinearMomentum
from Grid import GridHandlerGMSH
import Utils as utils
import dolfin as do
import torch as to
import numpy as np
import os
import copy

class Simulator(object):
	"""
	This class is responsible to carry out the simulation according to the
	input_file.json. It loads the mesh, builds the constitutive model, builds
	the linear momenum equation, defines the solver, defines the weak formulation
	the run the simulation.

	Parameters
	----------
	input_file : dict
		Dictionary extracted from the JSON file.
	"""
	def __init__(self, input_file):
		self.input_file = input_file

		# This is input_file to be saved
		self.input_file_to_be_saved = copy.deepcopy(input_file)

		# Output folder
		self.output_folder = input_file["output"]["path"]

		# Transient settings
		self.time_list = input_file["time_settings"]["time_list"]
		theta = input_file["time_settings"]["theta"]
	
		# Create mesh
		self.grid = GridHandlerGMSH(input_file["grid"]["name"], input_file["grid"]["path"])

		# Get number of elements
		self.n_elems = self.grid.mesh.num_cells()

		# Define constitutive model
		self.m = ConstitutiveModel(self.n_elems, input_file["constitutive_model"])

		# Build linear momentum equation handler
		self.eq = LinearMomentum(self.m, theta)

		# Define solver
		self.solver = self.define_solver()

		# Save input file
		filename = os.path.join(os.path.join(input_file["output"]["path"], "input_file.json"))
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		utils.save_json(self.input_file_to_be_saved, filename)

		# Define salt specific weight
		self.gravity = input_file["body_force"]["gravity"]
		self.density = input_file["body_force"]["density"]
		self.direction = input_file["body_force"]["direction"]

		# Create function spaces
		self.CG_3x1 = do.VectorFunctionSpace(self.grid.mesh, "CG", 1)
		self.DG_1x1 = do.FunctionSpace(self.grid.mesh, "DG", 0)
		self.DG_3x3 = do.TensorFunctionSpace(self.grid.mesh, "DG", 0)
		self.DG_6x6 = do.TensorFunctionSpace(self.grid.mesh, "DG", 0, shape=(6, 6))

		# Create tensor fields
		self.C0 = do.Function(self.DG_6x6)
		self.C1 = do.Function(self.DG_6x6)
		self.CT = do.Function(self.DG_6x6)
		self.eps_tot = do.Function(self.DG_3x3)
		self.eps_rhs = do.Function(self.DG_3x3)
		self.sigma = do.Function(self.DG_3x3)
		self.sigma_0 = do.Function(self.DG_3x3)

		# Define variational problem
		self.du = do.TrialFunction(self.CG_3x1)
		self.v = do.TestFunction(self.CG_3x1)
		self.ds = do.Measure("ds", domain=self.grid.mesh, subdomain_data=self.grid.get_boundaries())
		self.dx = do.Measure("dx", domain=self.grid.mesh, subdomain_data=self.grid.get_subdomains())
		self.normal = do.dot(self.v, do.FacetNormal(self.grid.mesh))
		self.n = do.FacetNormal(self.grid.mesh)

		# Create displacement vector
		self.u = do.Function(self.CG_3x1)
		self.u.rename("Displacement", "m")

		# Build body forces on RHS vector
		f_form = [0, 0, 0]
		f_form[self.direction] = self.gravity*self.density
		f_form = tuple(f_form)
		f = do.Constant(f_form)
		self.b_body = do.dot(f, self.v)*self.dx

	def run(self):
		"""
		Runs transient simulation.
		"""
		if self.input_file["simulation_settings"]["equilibrium"]["active"] == True:
			self.run_equilibrium()
			self.run_operation()
		else:
			self.run_simulation()

	def run_simulation(self):
		"""
		Runs simulation **without** solving the equilibrium condition.
		"""
		# Pseudo time
		t = self.time_list[0]
		t_final = self.time_list[-1]

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["operation"]["dt_max"]

		# Read number of time steps to skip before saving results
		n_skip = self.input_file["simulation_settings"]["operation"]["n_skip"]

		# Output folder
		operation_output_folder = os.path.join(self.input_file["output"]["path"], "operation")

		# Define Dirichlet boundary conditions
		bcs = self.define_dirichlet_bc(t)

		# Apply Neumann boundary conditions
		b_outer = self.apply_neumann_bc(t)

		# Build RHS vector
		b = do.assemble(self.b_body + b_outer)

		# Compute initial hardening
		for elem in self.eq.elems_ie:
			try:
				# Compute initial hardening parameter (alpha_0) based on initial stresses
				elem.compute_initial_hardening(self.eq.stress, Fvp_0=0.0)
				
				# Compute initial yield function values
				I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(*elem.extract_stress_components(self.eq.stress))
				_ = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)
				# Fvp.vector()[:] = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)

				print("Fvp: ", float(max(elem.Fvp)))
				print("alpha_min: ", float(min(elem.alpha)))
				print("alpha_max: ", float(max(elem.alpha)))
				print("alpha_avg: ", float(np.average(elem.alpha)))
				print()
			except:
				pass

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# Compute stress
		eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		self.eq.compute_stress_C0(eps_tot_torch)

		# Compute old ielastic strain rates
		self.eq.compute_eps_ie_rate()
		# m.compute_eps_ve_rate(0)

		# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.eq.update_eps_ie_rate_old()
		# m.update_eps_ve_rate_old()

		# Assign stress
		self.sigma.vector()[:] = self.eq.stress.flatten()

		# Create output file
		u_vtk = do.File(os.path.join(operation_output_folder, "vtk", "displacement", "displacement.pvd"))
		stress_vtk = do.File(os.path.join(operation_output_folder, "vtk", "stress", "stress.pvd"))

		# Save output fields
		u_vtk << (self.u, t)
		stress_vtk << (self.sigma, t)

		# Transient simulation
		n_step = 1
		while t < t_final:

			# Increment time
			t += dt

			# Define Dirichlet boundary conditions
			bcs = self.define_dirichlet_bc(t)

			# Apply Neumann boundary conditions
			b_outer = self.apply_neumann_bc(t)

			# Compute GT and BT matrix fields for viscoelastic elements
			self.eq.compute_GT_BT_ve(dt)

			# Iterative loop settings
			tol = 1e-7
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

				# Update stress of previous iteration (stress_k <-- stress)
				self.eq.update_stress()

				# Compute GT and BT matrix fields for inelastic elements
				self.eq.compute_GT_BT_ie(dt)

				# Compute CT
				self.eq.compute_CT(dt)
				self.CT.vector()[:] = to.flatten(self.eq.CT)

				# Compute right-hand side
				self.eq.compute_eps_rhs(dt)

				# Assign eps_rhs
				self.eps_rhs.vector()[:] = self.eq.eps_rhs.flatten()

				# Build rhs
				b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
				b = assemble(self.b_body + b_outer + b_rhs)

				# Build lhs
				a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
				A = assemble(a_form)

				# Solve linear system
				[bc.apply(A, b) for bc in bcs]
				self.solver.solve(A, self.u.vector(), b)
				# try:
				# except:
				# 	print("***** KABUM *****")
				# 	t = 2*t_final
				# 	break

				# Compute total strain
				self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
				eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

				# Compute stress
				self.eq.compute_stress(eps_tot_torch, dt)

				# Compute strain rates
				self.eq.compute_eps_ie_rate()
				self.eq.compute_eps_ve_rate(dt)

				# Compute error
				if self.eq.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k)
					eps_tot_flat = self.eps_tot.vector()[:]
					error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

				# Increment iteration counter
				ite += 1

			# Compute strains
			self.eq.compute_eps_ie(dt)
			self.eq.compute_eps_ve(dt)

			# Update old non-elastic strains
			self.eq.update_eps_ie_old()
			self.eq.update_eps_ve_old()

			# Update old non-elastic strain rates
			self.eq.update_eps_ie_rate_old()
			self.eq.update_eps_ve_rate_old()

			# Save displacement field
			u_vtk << (self.u, t)
			self.sigma.vector()[:] = self.eq.stress.flatten()
			stress_vtk << (self.sigma, t)

			# Save displacement field
			if n_step % n_skip == 0:
				u_vtk << (self.u, t)
				self.sigma.vector()[:] = self.eq.stress.flatten()
				stress_vtk << (self.sigma, t)

			# Print stuff
			print(n_step, f"{t_final/utils.hour}", t/utils.hour, ite, error)
			print()
			n_step += 1

	def run_operation(self):
		"""
		Runs transient simulation **after** the equilibrium condition. 
		"""
		# Pseudo time
		t = self.time_list[0]
		t_final = self.time_list[-1]

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["operation"]["dt_max"]

		# Read number of time steps to skip before saving results
		n_skip = self.input_file["simulation_settings"]["operation"]["n_skip"]

		# Output folder
		operation_output_folder = os.path.join(self.input_file["output"]["path"], "operation")

		# Define Dirichlet boundary conditions
		bcs = self.define_dirichlet_bc(t)

		# Apply Neumann boundary conditions
		b_outer = self.apply_neumann_bc(t)

		# Build RHS vector
		b = do.assemble(self.b_body + b_outer)

		# Compute initial hardening
		for elem in self.eq.m.elems_ie:
			try:
				# Compute initial hardening parameter (alpha_0) based on initial stresses
				elem.compute_initial_hardening(self.eq.stress, Fvp_0=0.0)
				
				# Compute initial yield function values
				I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(*elem.extract_stress_components(self.eq.stress))
				_ = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)
				# Fvp.vector()[:] = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)

				print("Fvp: ", float(max(elem.Fvp)))
				print("alpha_min: ", float(min(elem.alpha)))
				print("alpha_max: ", float(max(elem.alpha)))
				print("alpha_avg: ", float(np.average(elem.alpha)))
				print()
			except:
				pass

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# Compute stress
		eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		self.eq.compute_stress_C0(eps_tot_torch)

		# Compute old ielastic strain rates
		self.eq.compute_eps_ie_rate()
		# m.compute_eps_ve_rate(0)

		# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.eq.update_eps_ie_rate_old()
		# m.update_eps_ve_rate_old()

		# Assign stress
		self.sigma.vector()[:] = self.eq.stress.flatten()

		# Create output file
		u_vtk = do.File(os.path.join(operation_output_folder, "vtk", "displacement", "displacement.pvd"))
		stress_vtk = do.File(os.path.join(operation_output_folder, "vtk", "stress", "stress.pvd"))

		# Save output fields
		u_vtk << (self.u, t)
		stress_vtk << (self.sigma, t)

		# Transient simulation
		n_step = 1
		while t < t_final:

			# Increment time
			t += dt

			# Define Dirichlet boundary conditions
			bcs = self.define_dirichlet_bc(t)

			# Apply Neumann boundary conditions
			b_outer = self.apply_neumann_bc(t)

			# Compute GT and BT matrix fields for viscoelastic elements
			self.eq.compute_GT_BT_ve(dt)

			# Iterative loop settings
			tol = 1e-7
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

				# Update stress of previous iteration (stress_k <-- stress)
				self.eq.update_stress()

				# Compute GT and BT matrix fields for inelastic elements
				self.eq.compute_GT_BT_ie(dt)

				# Compute CT
				self.eq.compute_CT(dt)
				self.CT.vector()[:] = to.flatten(self.eq.CT)

				# Compute right-hand side
				self.eq.compute_eps_rhs(dt)

				# Assign eps_rhs
				self.eps_rhs.vector()[:] = self.eq.eps_rhs.flatten()

				# Build rhs
				b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
				b = do.assemble(self.b_body + b_outer + b_rhs)

				# Build lhs
				a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
				A = do.assemble(a_form)

				# Solve linear system
				[bc.apply(A, b) for bc in bcs]
				self.solver.solve(A, self.u.vector(), b)
				# try:
				# except:
				# 	print("***** KABUM *****")
				# 	t = 2*t_final
				# 	break

				# Compute total strain
				self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
				eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

				# Compute stress
				self.eq.compute_stress(eps_tot_torch, dt)

				# Compute strain rates
				self.eq.compute_eps_ie_rate()
				self.eq.compute_eps_ve_rate(dt)

				# Compute error
				if self.eq.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k)
					eps_tot_flat = self.eps_tot.vector()[:]
					error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

				# Increment iteration counter
				ite += 1

			# Compute strains
			self.eq.compute_eps_ie(dt)
			self.eq.compute_eps_ve(dt)

			# Update old non-elastic strains
			self.eq.update_eps_ie_old()
			self.eq.update_eps_ve_old()

			# Update old non-elastic strain rates
			self.eq.update_eps_ie_rate_old()
			self.eq.update_eps_ve_rate_old()

			# Save displacement field
			u_vtk << (self.u, t)
			self.sigma.vector()[:] = self.eq.stress.flatten()
			stress_vtk << (self.sigma, t)

			# Save displacement field
			if n_step % n_skip == 0:
				u_vtk << (self.u, t)
				self.sigma.vector()[:] = self.eq.stress.flatten()
				stress_vtk << (self.sigma, t)

			# Print stuff
			print(n_step, f"{t_final/utils.hour}", t/utils.hour, ite, error)
			print()
			n_step += 1

	def run_equilibrium(self):
		"""
		Runs simulation to determine the equilibrium condition. The equilibrium condition
		only considers the **elastic** and **viscoelastic** (if present) elements of the 
		constitutive model.
		"""
		# Pseudo time
		t_0 = self.time_list[0]
		t = t_0

		# Get maximum time step size
		dt = self.input_file["simulation_settings"]["equilibrium"]["dt_max"]

		# Output folder
		equilibrium_output_folder = os.path.join(self.input_file["output"]["path"], "equilibrium")

		# Define Dirichlet boundary conditions
		bcs = self.define_dirichlet_bc(t_0)
		
		# Apply Neumann boundary conditions
		b_outer = self.apply_neumann_bc(t_0)
		
		# Build RHS vector
		b = do.assemble(self.b_body + b_outer)

		# Initialize elastic stiffness matrix, C0
		self.C0.vector()[:] = to.flatten(self.eq.m.C0)

		# Build stiffness matrix
		a_form = do.inner(utils.dotdot(self.C0, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
		A = do.assemble(a_form)

		# Solve linear system
		[bc.apply(A, b) for bc in bcs]
		self.solver.solve(A, self.u.vector(), b)

		# Compute total strain
		self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))

		# Compute stress
		eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))
		self.eq.compute_stress_C0(eps_tot_torch)

		# Compute old viscoelastic strain rates
		self.eq.compute_eps_ve_rate(0)

		# Update viscoelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
		self.eq.update_eps_ve_rate_old()

		# Create output file
		u_vtk = do.File(os.path.join(equilibrium_output_folder, "vtk", "displacement", "displacement.pvd"))
		stress_vtk = do.File(os.path.join(equilibrium_output_folder, "vtk", "stress", "stress.pvd"))

		# Save output fields
		u_vtk << (self.u, t)
		stress_vtk << (self.sigma, t)

		# Initialize pseudo time control settings
		n_step = 1
		tol_time = self.input_file["simulation_settings"]["equilibrium"]["time_tol"]
		error_time = 2*tol_time
		eps_tot_old = utils.numpy2torch(self.eps_tot.vector()[:])

		while error_time > tol_time or n_step <= 2:

			# Increment time
			t += dt

			# Compute GT matrix field for viscoelastic elements
			self.eq.compute_GT_BT_ve(dt)

			# Iterative loop settings
			tol = 1e-7
			error = 2*tol
			ite = 0
			maxiter = 40

			while error > tol and ite < maxiter:

				# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
				eps_tot_k = utils.numpy2torch(self.eps_tot.vector()[:])

				# Update stress of previous iteration (stress_k <-- stress)
				self.eq.update_stress()

				# Compute CT
				self.eq.compute_CT(dt)
				self.CT.vector()[:] = to.flatten(self.eq.CT)

				# Compute right-hand side
				self.eq.compute_eps_rhs(dt)

				# Assign eps_rhs
				self.eps_rhs.vector()[:] = self.eq.eps_rhs.flatten()

				# Build rhs
				b_rhs = do.inner(utils.dotdot(self.CT, self.eps_rhs), utils.epsilon(self.v))*self.dx
				b = do.assemble(self.b_body + b_outer + b_rhs)

				# Build lhs
				a_form = do.inner(utils.dotdot(self.CT, utils.epsilon(self.du)), utils.epsilon(self.v))*self.dx
				A = do.assemble(a_form)

				# Solve linear system
				[bc.apply(A, b) for bc in bcs]
				self.solver.solve(A, self.u.vector(), b)

				# Compute total strain
				self.eps_tot.assign(utils.local_projection(utils.epsilon(self.u), self.DG_3x3))
				eps_tot_torch = utils.numpy2torch(self.eps_tot.vector()[:].reshape((self.n_elems, 3, 3)))

				# Compute stress
				self.eq.compute_stress(eps_tot_torch, dt)

				# Compute strain rates
				self.eq.compute_eps_ve_rate(dt)

				# Compute error
				if self.eq.theta == 1.0:
					error = 0.0
				else:
					eps_tot_k_flat = to.flatten(eps_tot_k)
					eps_tot_flat = self.eps_tot.vector()[:]
					error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

				# Increment iteration counter
				ite += 1

			# Compute strains
			self.eq.compute_eps_ve(dt)

			# Update old viscoelastic strains
			self.eq.update_eps_ve_old()

			# Update old viscoelastic strain rates
			self.eq.update_eps_ve_rate_old()

			# Compute time error (to check if steady state is achieved)
			eps_tot_flat = self.eps_tot.vector()[:]
			error_time = np.linalg.norm(eps_tot_old - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
			eps_tot_old = utils.numpy2torch(self.eps_tot.vector()[:])

			# Save displacement field
			u_vtk << (self.u, t)
			self.sigma.vector()[:] = self.eq.stress.flatten()
			stress_vtk << (self.sigma, t)

			# Print stuff
			print(n_step, t/utils.hour, ite, error, error_time)
			# print()
			n_step += 1



	def define_dirichlet_bc(self, t):
		"""
		Defines a list of Dirichlet boundary conditions to be applied to the
		linear system.

		Parameters
		----------
		t : float
			Time level of the simulation.

		Returns
		-------
		bcs : list[dolfin.fem.dirichletbc.DirichletBC]
			List containing the Dirichlet boundary conditions.
		"""
		i = 0
		bcs = []
		bc_dirichlet_list = []
		for boundary in self.input_file["boundary_conditions"]:
			if self.input_file["boundary_conditions"][boundary]["type"] == "dirichlet":
				bc_dirichlet_list.append(do.Expression("value", value=0, degree=1))
				component = self.input_file["boundary_conditions"][boundary]["component"]
				values = self.input_file["boundary_conditions"][boundary]["values"]
				bc_dirichlet_list[i].value = np.interp(t, self.time_list, values)
				bcs.append(do.DirichletBC(self.CG_3x1.sub(component), bc_dirichlet_list[i], self.grid.get_boundaries(), self.grid.get_boundary_tags(boundary)))
				i += 1
		return bcs

	def apply_neumann_bc(self, t):
		"""
		It reads all Neumann boundary conditions (external loads) applied to the geometry and
		builds the right-hand side vector.

		Parameters
		----------
		t : float
			Time level.

		Returns
		-------
		b_outer : ufl.form.Form
			This is the right-hand side vector corresponding to the Neumann boundary conditions,
			that is, it contains all the external loads applied to the geometry boundaries.
		"""
		i = 0
		bc_neumann_list = []
		for boundary in self.input_file["boundary_conditions"]:
			if self.input_file["boundary_conditions"][boundary]["type"] == "neumann":
				bc_direction = self.input_file["boundary_conditions"][boundary]["direction"]
				bc_density = self.input_file["boundary_conditions"][boundary]["density"]
				ref_position = self.input_file["boundary_conditions"][boundary]["reference_position"]
				bc_neumann_list.append(do.Expression(f"load_ref + rho*g*(H - x[{bc_direction}])", load_ref=0, rho=bc_density, g=self.gravity, H=ref_position, degree=1))
				values = self.input_file["boundary_conditions"][boundary]["values"]
				bc_neumann_list[i].load_ref = -np.interp(t, self.time_list, values)
				if i == 0: 	b_outer = bc_neumann_list[i]*self.normal*self.ds(self.grid.get_boundary_tags(boundary))
				else: 		b_outer += bc_neumann_list[i]*self.normal*self.ds(self.grid.get_boundary_tags(boundary))
				i += 1
		return b_outer


	def define_solver(self):
		"""
		Defines the solver for solving the linear system according to the specifications
		in the input_file.json.

		Returns
		-------
		solver : dolfin.cpp.la.KrylovSolver or dolfin.cpp.la.LUSolver
			The solver can be either an iterative solver (KrylovSolver) or a direct solver
			(LUSolver).
		"""
		solver_type = self.input_file["solver_settings"]["type"]
		if solver_type == "KrylovSolver":
			solver = do.KrylovSolver(
			                      	method = self.input_file["solver_settings"]["method"],
			                      	preconditioner = self.input_file["solver_settings"]["preconditioner"]
			                      )
			solver.parameters["relative_tolerance"] = self.input_file["solver_settings"]["relative_tolerance"]
		elif solver_type == "LU":
			solver = do.LUSolver(self.input_file["solver_settings"]["method"])
			solver.parameters["symmetric"] = self.input_file["solver_settings"]["symmetric"]
		else:
			raise Exception(f"Solver type {solver_type} not supported. Choose between KrylovSolver and LU.")
		return solver