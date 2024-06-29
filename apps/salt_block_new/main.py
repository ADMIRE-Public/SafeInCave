import os
import sys
sys.path.append(os.path.join("..", "..", "libs"))
from ConstitutiveModel import *
from Grid import GridHandlerGMSH
from Utils import *
from dolfin import *
import numpy as np
import torch as to
import pandas as pd
import copy
import json
import time

def compute_dt(Fvp_max):
	dt_max = 0.25*hour
	try:
		k3 = 1.4
		dt_min = 0.1*hour
		dt_star = dt_max/1.5

		k1 = (dt_max/(dt_max - dt_star)) - 1
		k2 = dt_max/(dt_max - dt_min)
		dt = dt_max - dt_max / (k2 + k1*to.exp(-k3*Fvp_max))
	except:
		dt = dt_max
		Fvp_max = -100
	print(f"dt = {dt/hour} h")
	print(f"Fvp: {float(Fvp_max)}")
	return float(dt)

def main():
	start_0 = time.time()

	# Read input file
	input_file = read_json("input_file_0.json")

	# This is input_file to be saved
	input_file_to_be_saved = copy.deepcopy(input_file)

	# Output folder
	output_folder = input_file["output"]["path"]
	
	# Create mesh
	g = GridHandlerGMSH(input_file["grid"]["name"], input_file["grid"]["path"])
	n_elems = g.mesh.num_cells()
	coordinates = g.mesh.coordinates()
	print("n_elems: ", n_elems)

	# Transient settings
	time_list = input_file["time_settings"]["time_list"]
	theta = input_file["time_settings"]["theta"]
	t = time_list[0]

	# Define salt specific weight
	gravity = input_file["body_force"]["gravity"]
	density = input_file["body_force"]["density"]
	direction = input_file["body_force"]["direction"]

	# Create function spaces
	CG_3x1 = VectorFunctionSpace(g.mesh, "CG", 1)
	DG_1x1 = FunctionSpace(g.mesh, "DG", 0)
	DG_3x3 = TensorFunctionSpace(g.mesh, "DG", 0)
	DG_6x6 = TensorFunctionSpace(g.mesh, "DG", 0, shape=(6, 6))

	# Create tensor fields
	C0 = Function(DG_6x6)
	C1 = Function(DG_6x6)
	CT = Function(DG_6x6)
	eps_tot = Function(DG_3x3)
	eps_rhs = Function(DG_3x3)
	sigma = Function(DG_3x3)
	sigma_0 = Function(DG_3x3)
	alpha_0 = Function(DG_1x1)
	alpha = Function(DG_1x1)
	Fvp = Function(DG_1x1)

	alpha.rename("Hardening parameter", "-")
	Fvp.rename("Yield function", "-")

	# Define variational problem
	du = TrialFunction(CG_3x1)
	v = TestFunction(CG_3x1)
	ds = Measure("ds", domain=g.mesh, subdomain_data=g.get_boundaries())
	dx = Measure("dx", domain=g.mesh, subdomain_data=g.get_subdomains())
	normal = dot(v, FacetNormal(g.mesh))
	n = FacetNormal(g.mesh)

	# Create displacement vector
	u = Function(CG_3x1)
	u.rename("Displacement", "m")

	# Define constitutive model
	m = ConstitutiveModelHandler(g, input_file)

	# # Add elastic element(s)
	# elems_e = get_list_of_elements_new(input_file["material_properties"], n_elems, element_class="Elastic")
	# for elem_e in elems_e:
	# 	m.add_elastic_element(elem_e)

	# # Add viscoelastic element
	# elems_ve = get_list_of_elements_new(input_file["material_properties"], n_elems, element_class="Viscoelastic")
	# for elem_ve in elems_ve:
	# 	m.add_viscoelastic_element(elem_ve)

	# # Add viscoelastic element
	# elems_ie = get_list_of_elements_new(input_file["material_properties"], n_elems, element_class="Inelastic")
	# for elem_ie in elems_ie:
	# 	m.add_inelastic_element(elem_ie)

	# # Initialize constitutive model
	# m.initialize()

	# Define Dirichlet boundary conditions
	i = 0
	bcs = []
	bc_dirichlet_list = []
	for boundary in input_file["boundary_conditions"]:
		if input_file["boundary_conditions"][boundary]["type"] == "dirichlet":
			bc_dirichlet_list.append(Expression("value", value=0, degree=1))
			component = input_file["boundary_conditions"][boundary]["component"]
			values = input_file["boundary_conditions"][boundary]["values"]
			bc_dirichlet_list[i].value = np.interp(t, time_list, values)
			bcs.append(DirichletBC(CG_3x1.sub(component), bc_dirichlet_list[i], g.get_boundaries(), g.get_boundary_tags(boundary)))
			i += 1

	# Apply Neumann boundary conditions
	i = 0
	bc_neumann_list = []
	for boundary in input_file["boundary_conditions"]:
		if input_file["boundary_conditions"][boundary]["type"] == "neumann":
			bc_direction = input_file["boundary_conditions"][boundary]["direction"]
			bc_density = input_file["boundary_conditions"][boundary]["density"]
			ref_position = input_file["boundary_conditions"][boundary]["reference_position"]
			bc_neumann_list.append(Expression(f"load_ref + rho*g*(H - x[{bc_direction}])", load_ref=0, rho=bc_density, g=gravity, H=ref_position, degree=1))
			values = input_file["boundary_conditions"][boundary]["values"]
			bc_neumann_list[i].load_ref = -np.interp(t, time_list, values)
			if i == 0: 	b_outer = bc_neumann_list[i]*normal*ds(g.get_boundary_tags(boundary))
			else: 		b_outer += bc_neumann_list[i]*normal*ds(g.get_boundary_tags(boundary))
			i += 1

	# Define solver
	if input_file["solver_settings"]["type"] == "KrylovSolver":
		solver = KrylovSolver(
		                      	method = input_file["solver_settings"]["method"],
		                      	preconditioner = input_file["solver_settings"]["preconditioner"]
		                      )
		solver.parameters["relative_tolerance"] = input_file["solver_settings"]["relative_tolerance"]

	elif input_file["solver_settings"]["type"] == "LU":
		solver = LUSolver(input_file["solver_settings"]["method"])
		solver.parameters["symmetric"] = input_file["solver_settings"]["symmetric"]

	# Build RHS vector
	f_form = [0, 0, 0]
	f_form[direction] = gravity*density
	f_form = tuple(f_form)
	f = Constant(f_form)
	b_body = dot(f, v)*dx
	b = assemble(b_body + b_outer)

	# Initialize elastic stiffness matrix, C0
	C0.vector()[:] = to.flatten(m.C0)
	C_aux = to.zeros_like(to.flatten(m.C0))
	for elem in m.elems_ve:
		C_aux += to.flatten(elem.C1)
	C1.vector()[:] = C_aux

	# Build stiffness matrix
	a_form = inner(dotdot(C0+1*C1, epsilon(du)), epsilon(v))*dx
	A = assemble(a_form)



	# Solve linear system
	[bc.apply(A, b) for bc in bcs]
	solver.solve(A, u.vector(), b)
	# solve(A, u.vector(), b, "petsc")

	# Compute total strain
	eps_tot.assign(local_projection(epsilon(u), DG_3x3))

	# Compute stress
	eps_tot_torch = to_tensor(eps_tot.vector()[:].reshape((n_elems, 3, 3)))
	m.compute_stress_C0(eps_tot_torch)

	for elem in m.elems_ie:
		try:
			# Compute initial hardening parameter (alpha_0) based on initial stresses
			elem.compute_initial_hardening(m.stress, Fvp_0=0.0)
			alpha.vector()[:] = elem.alpha
			
			# Compute initial yield function values
			I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(*elem.extract_stress_components(m.stress))
			Fvp.vector()[:] = elem.compute_Fvp(elem.alpha, I1_star, J2, Sr)

			print("Fvp: ", float(max(Fvp.vector()[:])))
			print("alpha_min: ", float(min(alpha.vector()[:])))
			print("alpha_max: ", float(max(alpha.vector()[:])))
			print("alpha_avg: ", float(np.average(alpha.vector()[:])))
			print()
		except:
			pass

	print(alpha.vector()[:])



	# Compute old inelastic strain rates
	m.compute_eps_ie_rate()
	m.compute_eps_ve_rate(0)

	# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
	m.update_eps_ie_rate_old()
	m.update_eps_ve_rate_old()


	# Create displacement file
	u_vtk = File(os.path.join(output_folder, "vtk", "displacement", "displacement.pvd"))
	Fvp_vtk = File(os.path.join(output_folder, "vtk", "Fvp", "Fvp.pvd"))
	alpha_vtk = File(os.path.join(output_folder, "vtk", "alpha", "alpha.pvd"))
	stress_vtk = File(os.path.join(output_folder, "vtk", "stress", "stress.pvd"))

	# Save displacement field
	u_vtk << (u, t)
	Fvp_vtk << (Fvp, t)
	alpha_vtk << (alpha, t)
	stress_vtk << (sigma, t)

	n_step = 1
	t_final = time_list[-1]
	stress_old = m.stress.clone()
	try:
		Fvp_max = max(m.elems_ie[0].Fvp)
	except:
		Fvp_max = -100
	dt = compute_dt(Fvp_max)
	print()

	# Start transient simulation
	while t < t_final:

		# Increment time
		t += dt

		# Update Dirichlet boundary conditions
		i = 0
		bcs = []
		for boundary in input_file["boundary_conditions"]:
			if input_file["boundary_conditions"][boundary]["type"] == "dirichlet":
				component = input_file["boundary_conditions"][boundary]["component"]
				values = input_file["boundary_conditions"][boundary]["values"]
				bc_dirichlet_list[i].value = -np.interp(t, time_list, values)
				bcs.append(DirichletBC(CG_3x1.sub(component), bc_dirichlet_list[i], g.get_boundaries(), g.get_boundary_tags(boundary)))
				i += 1

		# Update Neumann boundary conditions
		i = 0
		for boundary in input_file["boundary_conditions"]:
			if input_file["boundary_conditions"][boundary]["type"] == "neumann":
				values = input_file["boundary_conditions"][boundary]["values"]
				bc_neumann_list[i].load_ref = -np.interp(t, time_list, values)
				if i == 0: 	b_outer = bc_neumann_list[i]*normal*ds(g.get_boundary_tags(boundary))
				else: 		b_outer += bc_neumann_list[i]*normal*ds(g.get_boundary_tags(boundary))
				i += 1

		# Iterative loop settings
		tol = 1e-7
		error = 2*tol
		ite = 0
		maxiter = 40

		while error > tol and ite < maxiter:

			# Update total strain of previous iteration (eps_tot_k <-- eps_tot)
			eps_tot_k = to_tensor(eps_tot.vector()[:])

			# Update stress of previous iteration (stress_k <-- stress)
			m.update_stress()

			# Compute GT matrix field
			m.compute_GT_BT_ve(dt)
			m.compute_GT_BT_ie(dt)

			# Compute CT
			m.compute_CT(dt)
			CT.vector()[:] = to.flatten(m.CT)

			# Compute right-hand side
			m.compute_eps_rhs(dt)

			# Assign eps_rhs
			eps_rhs.vector()[:] = m.eps_rhs.flatten()

			# Build rhs
			b_rhs = inner(dotdot(CT, eps_rhs), epsilon(v))*dx
			b = assemble(b_body + b_outer + b_rhs)

			# Build lhs
			a_form = inner(dotdot(CT, epsilon(du)), epsilon(v))*dx
			A = assemble(a_form)

			# Solve linear system
			[bc.apply(A, b) for bc in bcs]
			solver.solve(A, u.vector(), b)

			# Compute total strain
			eps_tot.assign(local_projection(epsilon(u), DG_3x3))
			eps_tot_torch = to_tensor(eps_tot.vector()[:].reshape((n_elems, 3, 3)))

			# Compute stress
			m.compute_stress(eps_tot_torch, dt)

			# Increment internal variables of non-elastic elements
			m.increment_internal_variables(dt)

			# Compute strain rates
			m.compute_eps_ie_rate()
			m.compute_eps_ve_rate(dt)

			# Compute error
			if theta == 1.0:
				error = 0.0
			else:
				eps_tot_k_flat = to.flatten(eps_tot_k)
				eps_tot_flat = eps_tot.vector()[:]
				error = np.linalg.norm(eps_tot_k_flat - eps_tot_flat) / np.linalg.norm(eps_tot_flat)

			# Increment iteration counter
			ite += 1

		# Update qsi_old
		m.update_internal_variables()

		# Compute strains
		m.compute_eps_ie(dt)
		m.compute_eps_ve(dt)

		# Update old non-elastic strains
		m.update_eps_ie_old()
		m.update_eps_ve_old()

		# Update old non-elastic strain rates
		m.update_eps_ie_rate_old()
		m.update_eps_ve_rate_old()

		for elem in m.elems_ie:
			try:
				Fvp.vector()[:] = elem.Fvp
				alpha.vector()[:] = elem.alpha
				Fvp_vtk << (Fvp, t)
				alpha_vtk << (alpha, t)
			except:
				pass


		# Compute time step size
		try:
			Fvp_max = max(m.elems_ie[0].Fvp)
		except:
			Fvp_max = -100
		dt = compute_dt(Fvp_max)

		# Update stress
		stress_old = m.stress.clone()

		# Save displacement field
		u_vtk << (u, t)
		sigma.vector()[:] = m.stress.flatten()
		stress_vtk << (sigma, t)

		# Print stuff
		print(n_step, f"{t_final/hour}", t/hour, ite, error)
		print()
		n_step += 1

	# Save inputs
	save_json(input_file_to_be_saved, os.path.join(output_folder, "input_file.json"))

	# Print simulation CPU time
	final = time.time()
	print(f"Time: {final - start_0} seconds.")

if __name__ == '__main__':
	main()