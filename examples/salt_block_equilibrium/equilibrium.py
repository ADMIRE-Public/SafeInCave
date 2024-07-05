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

def compute_equilibrium(g, m, solver, input_file):
	# Output folder
	output_folder = os.path.join(input_file["output"]["path"], "equilibrium")

	# Get number of elements
	n_elems = g.mesh.num_cells()
	
	# Transient settings
	time_list = input_file["time_settings"]["time_list"]
	theta = input_file["time_settings"]["theta"]
	t = time_list[0]

	# Get maximum time step size
	dt = input_file["simulation_settings"]["equilibrium"]["dt_max"]

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

	# Define Dirichlet boundary conditions
	i = 0
	bcs = []
	bc_dirichlet_list = []
	for boundary in input_file["boundary_conditions"]:
		if input_file["boundary_conditions"][boundary]["type"] == "dirichlet":
			bc_dirichlet_list.append(Expression("value", value=0, degree=1))
			component = input_file["boundary_conditions"][boundary]["component"]
			values = input_file["boundary_conditions"][boundary]["values"]
			bc_dirichlet_list[i].value = values[0]
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
			bc_neumann_list[i].load_ref = -values[0]
			if i == 0: 	b_outer = bc_neumann_list[i]*normal*ds(g.get_boundary_tags(boundary))
			else: 		b_outer += bc_neumann_list[i]*normal*ds(g.get_boundary_tags(boundary))
			i += 1


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
	a_form = inner(dotdot(C0, epsilon(du)), epsilon(v))*dx
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

	# Compute old viscoelastic strain rates
	m.compute_eps_ve_rate(0)

	# Update viscoelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
	m.update_eps_ve_rate_old()

	# Create output file
	u_vtk = File(os.path.join(output_folder, "vtk", "displacement", "displacement.pvd"))
	stress_vtk = File(os.path.join(output_folder, "vtk", "stress", "stress.pvd"))

	# Save output fields
	u_vtk << (u, t)
	stress_vtk << (sigma, t)

	n_step = 1
	tol_time = input_file["simulation_settings"]["equilibrium"]["time_tol"]
	error_time = 2*tol_time
	eps_tot_old = to_tensor(eps_tot.vector()[:])
	while error_time > tol_time or n_step <= 2:

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

			# Compute strain rates
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

		# Compute strains
		m.compute_eps_ve(dt)

		# Update old viscoelastic strains
		m.update_eps_ve_old()

		# Update old viscoelastic strain rates
		m.update_eps_ve_rate_old()

		# Compute time error (to check if steady state is achieved)
		eps_tot_flat = eps_tot.vector()[:]
		error_time = np.linalg.norm(eps_tot_old - eps_tot_flat) / np.linalg.norm(eps_tot_flat)
		eps_tot_old = to_tensor(eps_tot.vector()[:])

		# Save displacement field
		u_vtk << (u, t)
		sigma.vector()[:] = m.stress.flatten()
		stress_vtk << (sigma, t)

		# Print stuff
		print(n_step, t/hour, ite, error, error_time)
		print()
		n_step += 1


		if input_file["simulation_settings"]["equilibrium"]["active"] == False:
			return u

	return u

if __name__ == '__main__':
	main()