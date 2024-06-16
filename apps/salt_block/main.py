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
	dt_max = 0.5*hour
	try:
		k3 = 1.4
		dt_min = 0.01*hour
		dt_star = dt_max/1.5

		k1 = (dt_max/(dt_max - dt_star)) - 1
		k2 = dt_max/(dt_max - dt_min)
		dt = dt_max - dt_max / (k2 + k1*to.exp(-k3*Fvp_max))
		# print(float(Fvp_max), float(dt), dt_max)
	except:
		dt = dt_max
		Fvp_max = -100
	print(f"dt = {dt/hour} h")
	print(f"Fvp: {float(Fvp_max)}")
	return float(dt)

def main():
	start_0 = time.time()

	# Read input_bc
	# input_bc = read_json("input_bc_1.json")
	input_bc = read_json("input_bc_0.json")
	time_refined = np.array(input_bc["Time"]["timeList"])
	sigma_axial = np.array(input_bc["sigma_axial"])
	sigma_radial = np.array(input_bc["sigma_radial"])

	# time_list = np.linspace(time_refined[0], time_refined[-1], 400)
	time_list = time_refined

	def get_loads(t):
		s_axi = -np.interp(t, time_refined, sigma_axial)
		s_rad = -np.interp(t, time_refined, sigma_radial)
		return s_rad, s_axi
	
	# Read input model
	input_model = read_json("input_model.json")

	# This is input_model to be saved
	input_model_to_be_saved = copy.deepcopy(input_model)

	# Create mesh
	# g = GridHandlerGMSH("geom", os.path.join("..", "..", "grids", "sugar_cube_0"))
	g = GridHandlerGMSH("geom", os.path.join("..", "..", "grids", "cube_0"))
	n_elems = g.mesh.num_cells()
	coordinates = g.mesh.coordinates()
	print("n_elems: ", n_elems)

	# Create function spaces
	VS = VectorFunctionSpace(g.mesh, "CG", 1)
	TS = TensorFunctionSpace(g.mesh, "DG", 0)
	P0 = FunctionSpace(g.mesh, "DG", 0)
	V_DG_6x6 = TensorFunctionSpace(g.mesh, "DG", 0, shape=(6, 6))

	# Create tensor fields
	C0 = Function(V_DG_6x6)
	C1 = Function(V_DG_6x6)
	CT = Function(V_DG_6x6)
	eps_tot = Function(TS)
	eps_rhs = Function(TS)
	sigma = Function(TS)
	sigma_0 = Function(TS)
	alpha_0 = Function(P0)
	alpha = Function(P0)
	Fvp = Function(P0)

	alpha.rename("Hardening parameter", "-")
	Fvp.rename("Yield function", "-")

	# Define boundary condition
	bcs = []
	bcs.append(DirichletBC(VS.sub(0), Constant(0.0), g.get_boundaries(), g.get_boundary_tags("WEST")))
	bcs.append(DirichletBC(VS.sub(1), Constant(0.0), g.get_boundaries(), g.get_boundary_tags("SOUTH")))
	bcs.append(DirichletBC(VS.sub(2), Constant(0.0), g.get_boundaries(), g.get_boundary_tags("BOTTOM")))

	# Define variational problem
	du = TrialFunction(VS)
	v = TestFunction(VS)
	ds = Measure("ds", domain=g.mesh, subdomain_data=g.get_boundaries())
	dx = Measure("dx", domain=g.mesh, subdomain_data=g.get_subdomains())
	normal = dot(v, FacetNormal(g.mesh))
	n = FacetNormal(g.mesh)

	# Define boundary conditions
	

	# Create displacement vector
	u = Function(VS)
	u_0 = Function(VS)
	delta_u = Function(VS)
	delta_u.rename("Displacement", "m")


	# Transient settings
	t = time_list[0]
	theta = input_bc["Time"]["theta"]

	# Output folder
	output_folder = os.path.join("output", "case_IMP")
	print(output_folder)
	print()

	# Define constitutive model
	m = ConstitutiveModelHandler(theta, n_elems)

	# Add elastic element(s)
	elems_e = get_list_of_elements(input_model, n_elems, element_class="Elastic")
	for elem_e in elems_e:
		m.add_elastic_element(elem_e)

	# Add viscoelastic element
	elems_ve = get_list_of_elements(input_model, n_elems, element_class="Viscoelastic")
	for elem_ve in elems_ve:
		m.add_viscoelastic_element(elem_ve)

	# Add viscoelastic element
	elems_ie = get_list_of_elements(input_model, n_elems, element_class="Inelastic")
	for elem_ie in elems_ie:
		m.add_inelastic_element(elem_ie)

	
	# Initialize constitutive model
	m.initialize()

	# Initialize elastic stiffness matrix, C0
	C0.vector()[:] = to.flatten(m.C0)
	C_aux = to.zeros_like(to.flatten(m.C0))
	for elem in m.elems_ve:
		C_aux += to.flatten(elem.C1)
	C1.vector()[:] = C_aux

	# Build RHS vector
	f = Constant((0, 0, 0))
	b_body = dot(f, v)*dx
	s_rad, s_axi = get_loads(t)
	b_outer  = s_rad*normal*ds(g.get_boundary_tags("EAST"))
	b_outer += s_rad*normal*ds(g.get_boundary_tags("NORTH"))
	b_outer += s_axi*normal*ds(g.get_boundary_tags("TOP"))
	b = assemble(b_body + b_outer)

	s_rad, s_axi = get_loads(t)
	print("Axial load: %.2f MPa"%(s_axi/MPa))
	print("Radial load: %.2f MPa"%(s_rad/MPa))

	# Build stiffness matrix
	a_form = inner(dotdot(C0+1*C1, epsilon(du)), epsilon(v))*dx
	A = assemble(a_form)

	# Solve linear system
	[bc.apply(A, b) for bc in bcs]
	solver = KrylovSolver('cg', 'sor')
	# solver.parameters['absolute_tolerance'] = 1e-10
	solver.parameters['relative_tolerance'] = 1e-12
	solver.solve(A, u.vector(), b)
	# solve(A, u.vector(), b, "petsc")

	# Assign initial displacement field
	u_0.assign(u)

	# Compute total strain
	eps_tot.assign(local_projection(epsilon(u), TS))

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
	delta_u.vector()[:] = u.vector()[:] #- u_0.vector()[:]
	u_vtk << (delta_u, t)
	Fvp_vtk << (Fvp, t)
	alpha_vtk << (alpha, t)
	stress_vtk << (sigma, t)

	# for i in range(1, len(time_list)):
	# for i in range(1, 10):
	n_step = 1
	t_final = time_list[-1]
	# t_final = 1.5*hour
	stress_old = m.stress.clone()
	try:
		Fvp_max = max(m.elems_ie[0].Fvp)
	except:
		Fvp_max = -100
	dt = compute_dt(Fvp_max)
	print()
	while t < t_final:

		# dt = compute_dt(m.elems_ie[0].Fvp)
		t += dt

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
			s_rad, s_axi = get_loads(t)
			b_outer  = s_rad*normal*ds(g.get_boundary_tags("EAST"))
			b_outer += s_rad*normal*ds(g.get_boundary_tags("NORTH"))
			b_outer += s_axi*normal*ds(g.get_boundary_tags("TOP"))
			b = assemble(b_body + b_outer + b_rhs)

			# Build lhs
			a_form = inner(dotdot(CT, epsilon(du)), epsilon(v))*dx
			A = assemble(a_form)

			# Solve linear system
			[bc.apply(A, b) for bc in bcs]
			solver.solve(A, u.vector(), b)

			# Compute total strain
			eps_tot.assign(local_projection(epsilon(u), TS))
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

		s_rad, s_axi = get_loads(t)
		print("Axial load: %.2f MPa"%(s_axi/MPa))
		print("Radial load: %.2f MPa"%(s_rad/MPa))
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
		# dt = compute_dt_2(m.elems_ie[0].Fvp, m.stress, stress_old, dt)
		dt = compute_dt(Fvp_max)
		# dt = compute_dt(m.elems_ie[0].Fvp)

		# Update stress
		stress_old = m.stress.clone()

		# Save displacement field
		delta_u.vector()[:] = u.vector()[:] #- u_0.vector()[:]
		u_vtk << (delta_u, t)
		sigma.vector()[:] = m.stress.flatten()
		stress_vtk << (sigma, t)

		# Print stuff
		print(n_step, f"{t_final/hour}", t/hour, ite, error)
		print()
		n_step += 1

	# Save inputs
	save_json(input_bc, os.path.join(output_folder, "input_bc.json"))
	save_json(input_model_to_be_saved, os.path.join(output_folder, "input_model.json"))

	# Print simulation CPU time
	final = time.time()
	print(f"Time: {final - start_0} seconds.")

if __name__ == '__main__':
	main()