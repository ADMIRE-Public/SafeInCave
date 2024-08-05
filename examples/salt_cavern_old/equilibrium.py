import os
import sys
sys.path.append(os.path.join("..", "..", "libs"))
from Utils import *
from dolfin import *
import numpy as np
import torch as to

def compute_equilibrium(m, g, input_bc, output_folder):

	# Output folder
	output_folder = os.path.join(output_folder, "equilibrium")
	print(output_folder)
	print()

	t = 0
	dt = 2.0*hour

	time_refined = np.array(input_bc["Time"]["timeList"])
	gas_pressure = np.array(input_bc["gas_pressure"])

	def get_pressure(t):
		return gas_pressure[0]

	n_elems = g.mesh.num_cells()
	coordinates = g.mesh.coordinates()
	H = coordinates[:,2].max()
	print("n_elems: ", n_elems)

	# Create function spaces
	VS = VectorFunctionSpace(g.mesh, "CG", 1)
	TS = TensorFunctionSpace(g.mesh, "DG", 0)
	V_DG_6x6 = TensorFunctionSpace(g.mesh, "DG", 0, shape=(6, 6))

	# Create tensor fields
	C0 = Function(V_DG_6x6)
	C1 = Function(V_DG_6x6)
	CT = Function(V_DG_6x6)
	eps_tot = Function(TS)
	eps_rhs = Function(TS)
	sigma = Function(TS)

	# Define boundary condition
	bcs = []
	bcs.append(DirichletBC(VS.sub(0), Constant(0.0), g.get_boundaries(), g.get_boundary_tags("West")))
	bcs.append(DirichletBC(VS.sub(1), Constant(0.0), g.get_boundaries(), g.get_boundary_tags("South")))
	bcs.append(DirichletBC(VS.sub(2), Constant(0.0), g.get_boundaries(), g.get_boundary_tags("Bottom")))

	# Define variational problem
	du = TrialFunction(VS)
	v = TestFunction(VS)
	ds = Measure("ds", domain=g.mesh, subdomain_data=g.get_boundaries())
	dx = Measure("dx", domain=g.mesh, subdomain_data=g.get_subdomains())
	normal = dot(v, FacetNormal(g.mesh))
	n = FacetNormal(g.mesh)

	# Create displacement vector
	u = Function(VS)
	u.rename("Displacement", "m")

	# Initialize elastic stiffness matrix, C0
	C0.vector()[:] = to.flatten(m.C0)

	# Define salt specific weight
	rho_salt = 2000.0
	gravity = -1*9.8
	gamma_salt = 1*rho_salt*gravity

	# Define overburden and sideburden
	sigma_v = -10.0*MPa
	sigma_h = -10.0*MPa
	over_burden = Expression("sv", sv=sigma_v, degree=1)
	side_burden = Expression("sh + gamma*(H - x[2])", sh=sigma_h, gamma=gamma_salt, H=H, degree=1)
	# side_burden = Expression("sh + gamma*(H - x[2])", sh=sigma_h, gamma=0, H=H, degree=1)

	# Build RHS vector
	f = Constant((0, 0, gamma_salt))
	b_body = dot(f, v)*dx
	b_outer  = side_burden*normal*ds(g.get_boundary_tags("East"))
	b_outer += side_burden*normal*ds(g.get_boundary_tags("North"))
	b_outer += over_burden*normal*ds(g.get_boundary_tags("Top"))
	b_wall = get_pressure(t)*normal*ds(g.get_boundary_tags("Cavern"))
	b = assemble(b_body + b_outer + b_wall)

	print("Gas pressure: %.2f MPa"%(get_pressure(t)/MPa))

	# Build stiffness matrix
	a_form = inner(dotdot(C0, epsilon(du)), epsilon(v))*dx
	A = assemble(a_form)

	# Solve linear system
	[bc.apply(A, b) for bc in bcs]
	solver = KrylovSolver('cg', 'sor')
	# solver.parameters['absolute_tolerance'] = 1e-10
	solver.parameters['relative_tolerance'] = 1e-12
	solver.solve(A, u.vector(), b)
	# solve(A, u.vector(), b, "petsc")

	# Compute total strain
	eps_tot.assign(local_projection(epsilon(u), TS))

	# Compute stress
	eps_tot_torch = to_tensor(eps_tot.vector()[:].reshape((n_elems, 3, 3)))
	m.compute_stress_C0(eps_tot_torch)

	# Compute old inelastic strain rates
	m.compute_eps_ve_rate(0)

	# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
	m.update_eps_ve_rate_old()

	# Create displacement file
	u_vtk = File(os.path.join(output_folder, "u", "u.pvd"))
	stress_vtk = File(os.path.join(output_folder, "stress", "stress.pvd"))

	# Save displacement field
	u_vtk << (u, t)
	stress_vtk << (sigma, t)

	n_step = 1
	tol_time = 1e-4
	error_time = 2*tol_time
	eps_tot_old = to_tensor(eps_tot.vector()[:])
	while error_time > tol_time or n_step <= 2:

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

			# Compute CT
			m.compute_CT(dt)
			CT.vector()[:] = to.flatten(m.CT)


			# Compute right-hand side
			m.compute_eps_rhs(dt)

			# Assign eps_rhs
			eps_rhs.vector()[:] = m.eps_rhs.flatten()

			# Build rhs
			b_rhs = inner(dotdot(CT, eps_rhs), epsilon(v))*dx
			b_wall = get_pressure(t)*normal*ds(g.get_boundary_tags("Cavern"))
			b = assemble(b_body + b_outer + b_wall + b_rhs)

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

			# Compute strain rates
			m.compute_eps_ve_rate(dt)

			# Compute error
			if m.theta == 1.0:
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
		m.compute_eps_ve(dt)

		# Update old non-elastic strains
		m.update_eps_ve_old()

		# Update old non-elastic strain rates
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
		n_step += 1

	return u