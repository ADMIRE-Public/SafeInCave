import os
import sys
sys.path.append(os.path.join("..", "..", "libs"))
from Utils import *
from dolfin import *
import numpy as np
import torch as to

def compute_dt(Fvp_max):
	dt_max = 0.1*hour
	try:
		k3 = 0.5
		dt_min = 0.005*hour
		dt_star = dt_max/1.5

		k1 = (dt_max/(dt_max - dt_star)) - 1
		k2 = dt_max/(dt_max - dt_min)
		dt = dt_max - dt_max / (k2 + k1*to.exp(-k3*Fvp_max))
		# print(float(Fvp_max), float(dt), dt_max)
	except:
		dt = dt_max
		Fvp_max = -100
	# dt = 0.00001*hour
	print(f"dt = {dt/hour} h")
	print(f"Fvp: {float(Fvp_max)}")
	return float(dt)

def compute_q(stress_torch):
	s_xx = stress_torch[:,0,0]
	s_yy = stress_torch[:,1,1]
	s_zz = stress_torch[:,2,2]
	s_xy = stress_torch[:,0,1]
	s_xz = stress_torch[:,0,2]
	s_yz = stress_torch[:,1,2]
	q_vm = to.sqrt( 0.5*( (s_xx - s_yy)**2 + (s_xx - s_zz)**2 + (s_yy - s_zz)**2 + 6*(s_xy**2 + s_xz**2 + s_yz**2) ) )
	return q_vm



def run_simulation(m, g, u, input_bc, output_folder):

	# Output folder
	output_folder = os.path.join(output_folder, "operation")
	print(output_folder)
	print()

	time_refined = np.array(input_bc["Time"]["timeList"])
	gas_pressure = np.array(input_bc["gas_pressure"])

	t = 0
	t_final = time_refined[-1]

	def get_pressure(t):
		return np.interp(t, time_refined, gas_pressure)

	n_elems = g.mesh.num_cells()
	coordinates = g.mesh.coordinates()
	H = coordinates[:,2].max()
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
	alpha = Function(P0)
	Fvp = Function(P0)
	q_field = Function(P0)
	q_field.rename("Von Mises Stress", "(MPa)")

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

	# Define salt specific weight
	rho_salt = 2000.0
	gravity = -1*9.8
	gamma_salt = 1*rho_salt*gravity

	# Define overburden and sideburden
	sigma_v = -10.0*MPa
	sigma_h = -10.0*MPa
	over_burden = Expression("sv", sv=sigma_v, degree=1)
	side_burden = Expression("sh + gamma*(H - x[2])", sh=sigma_h, gamma=gamma_salt, H=H, degree=1)
	p_wall = Expression("p", p=0, degree=1)
	# side_burden = Expression("sh + gamma*(H - x[2])", sh=sigma_h, gamma=0, H=H, degree=1)

	# Build RHS vector
	f = Constant((0, 0, gamma_salt))
	b_body = dot(f, v)*dx
	b_outer  = side_burden*normal*ds(g.get_boundary_tags("East"))
	b_outer += side_burden*normal*ds(g.get_boundary_tags("North"))
	b_outer += over_burden*normal*ds(g.get_boundary_tags("Top"))

	print("Gas pressure: %.2f MPa"%(get_pressure(t)/MPa))

	# Define solver
	solver = KrylovSolver('cg', 'sor')
	solver.parameters['relative_tolerance'] = 1e-12

	# Compute total strain
	eps_tot.assign(local_projection(epsilon(u), TS))

	for elem in m.elems_ie:
		try:
			# Compute initial hardening parameter (alpha_0) based on initial stresses
			elem.compute_initial_hardening(m.stress, Fvp_0=0)
			alpha.vector()[:] = elem.alpha
			
			# Compute initial yield function values
			I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(*elem.extract_stress_components(m.stress))

			print("Fvp: ", float(max(Fvp.vector()[:])))
			print("alpha_min: ", float(min(alpha.vector()[:])))
			print("alpha_max: ", float(max(alpha.vector()[:])))
			print("alpha_avg: ", float(np.average(alpha.vector()[:])))
			print()
		except:
			pass

	print(alpha.vector()[:])
	# for i, value in enumerate(alpha.vector()[:]):
	# 	print(i, value)



	# Compute old inelastic strain rates
	m.compute_eps_ie_rate()
	# m.compute_eps_ve_rate(0)

	# Update inelastic strain rate (Warning! Do NOT update eps_ie here, because this is wrong!)
	m.update_eps_ie_rate_old()
	# m.update_eps_ve_rate_old()

	# Assign stress
	sigma.vector()[:] = m.stress.flatten()

	try:
		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = m.elems_ie[0].extract_stress_components(m.stress)
		I1, I2, I3, J2, J3, Sr, I1_star = m.elems_ie[0].compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
		I1_field.vector()[:] = I1_star
		J2_field.vector()[:] = J2
	except:
		pass


	# Create displacement file
	u_vtk = File(os.path.join(output_folder, "u", "u.pvd"))
	q_vtk = File(os.path.join(output_folder, "von_mises_stress", "q.pvd"))
	stress_vtk = File(os.path.join(output_folder, "stress", "stress.pvd"))

	# Save displacement field
	q_vm = compute_q(m.stress)
	q_field.vector()[:] = q_vm.flatten()
	u_vtk << (u, t)
	q_vtk << (q_field, t)
	stress_vtk << (sigma, t)

	n_step = 1
	stress_old = m.stress.clone()
	try:
		Fvp_max = max(m.elems_ie[0].Fvp)
	except:
		Fvp_max = -100
	dt = compute_dt(Fvp_max)
	print()
	# t_final = 5*hour
	while t < t_final:

		# dt = compute_dt(m.elems_ie[0].Fvp)
		t += dt

		p_wall.p = get_pressure(t)

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

			# for i, G in enumerate(m.elems_ie[0].G):
			# 	print(i, G)

			# Compute CT
			m.compute_CT(dt)
			CT.vector()[:] = to.flatten(m.CT)

			# Compute right-hand side
			m.compute_eps_rhs(dt)

			# Assign eps_rhs
			eps_rhs.vector()[:] = m.eps_rhs.flatten()

			# Build rhs
			b_rhs = inner(dotdot(CT, eps_rhs), epsilon(v))*dx
			b_wall = p_wall*normal*ds(g.get_boundary_tags("Cavern"))
			b = assemble(b_body + b_outer + b_wall + b_rhs)

			# for value in b.get_local():
			# 	print(value)

			# Build lhs
			a_form = inner(dotdot(CT, epsilon(du)), epsilon(v))*dx
			A = assemble(a_form)

			# Solve linear system
			[bc.apply(A, b) for bc in bcs]
			try:
				solver.solve(A, u.vector(), b)
			except:
				print("RUIM")
				t = 2*t_final
				break

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
		m.compute_eps_ie(dt)
		m.compute_eps_ve(dt)

		# Update old non-elastic strains
		m.update_eps_ie_old()
		m.update_eps_ve_old()

		# Update old non-elastic strain rates
		m.update_eps_ie_rate_old()
		m.update_eps_ve_rate_old()

		print("Gas pressure: %.2f MPa"%(get_pressure(t)/MPa))
		for elem in m.elems_ie:
			try:
				Fvp.vector()[:] = elem.Fvp
				alpha.vector()[:] = elem.alpha
				s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = elem.extract_stress_components(m.stress)
				I1, I2, I3, J2, J3, Sr, I1_star = elem.compute_stress_invariants(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
				# I1_field.vector()[:] = I1_star
				# J2_field.vector()[:] = J2
				# Fvp_vtk << (Fvp, t)
				# alpha_vtk << (alpha, t)
				# I1_vtk << (I1_field, t)
				# J2_vtk << (J2_field, t)
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

		# Compute von Mises stress
		q_vm = compute_q(m.stress)
		q_field.vector()[:] = q_vm.flatten()

		# Save displacement field
		if n_step%10 == 0:
		# if n_step%50 == 0:
			u_vtk << (u, t)
			sigma.vector()[:] = m.stress.flatten()
			stress_vtk << (sigma, t)
			q_vtk << (q_field, t)

		# Print stuff
		print(n_step, f"{t_final/hour}", t/hour, ite, error)
		print()
		n_step += 1