import safeincave as sf
import safeincave.Utils as ut
import torch as to
import matplotlib.pyplot as plt

MPa = 1e6

def main():
	n_elems = 1
	mat = sf.Material(n_elems)

	# Choose time integration method
	theta = 0.0

	# Elastic element
	E = 102e9*to.ones(n_elems)
	nu = 0.3*to.ones(n_elems)
	spring_0 = sf.Spring(E, nu, "spring")

	# Viscoelastic element
	eta = 105e11*to.ones(n_elems)
	E = 10e9*to.ones(n_elems)
	nu = 0.32*to.ones(n_elems)
	kelvin = sf.Viscoelastic(eta, E, nu, "kelvin")

	# Create creep
	A = 1.9e-20*to.ones(n_elems)
	Q = 51600*to.ones(n_elems)
	n = 3.0*to.ones(n_elems)
	creep_ds = sf.DislocationCreep(A, Q, n, "ds_creep")

	# Create pressure solution creep
	A = 1.29e-19*to.ones(n_elems)
	Q = 13184*to.ones(n_elems)
	d = 0.01*to.ones(n_elems)
	creep_ps = sf.PressureSolutionCreep(A, d, Q, "ps_creep")

	# Create Desai's viscoplastic model
	mu_1 = 5.3665857009859815e-11*to.ones(n_elems)
	N_1 = 3.1*to.ones(n_elems)
	n = 3.0*to.ones(n_elems)
	a_1 = 1.965018496922832e-05*to.ones(n_elems)
	eta = 0.8275682807874163*to.ones(n_elems)
	beta_1 = 0.0048*to.ones(n_elems)
	beta = 0.995*to.ones(n_elems)
	m = -0.5*to.ones(n_elems)
	gamma = 0.095*to.ones(n_elems)
	alpha_0 = 0.0022*to.ones(n_elems)
	sigma_t = 5.0*to.ones(n_elems)
	desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")

	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_ds)
	mat.add_to_non_elastic(creep_ps)
	# mat.add_to_non_elastic(desai)

	for elem_e in mat.elems_e:
		elem_e.initialize()

	# Define temperature
	Temp = 298
	T = Temp*to.ones(n_elems)

	# Define time_list
	t_0 = 0
	t_f = 3*ut.hour
	n_steps = 100
	time_list = to.linspace(t_0, t_f, n_steps)

	# Allocate eps_tot_list
	eps_tot_list = to.zeros((n_steps, n_elems, 3, 3))

	# Initialize stress
	stress_hist = to.zeros((n_steps, n_elems, 3, 3))
	stress_hist[:,:,0,0] = -5*MPa
	stress_hist[:,:,1,1] = -5*MPa
	stress_hist[:,:,2,2] = -15.5*MPa


	# Loop over time
	for i in range(1, n_steps):
		stress = stress_hist[i]
		dt = time_list[i] - time_list[i-1]
		phi1 = dt*theta
		phi2 = dt*(1 - theta)

		# Compute elastic strains
		for elem_e in mat.elems_e:
			elem_e.compute_eps_e(stress)

		# Initialize non-elastic strains of iteration k
		for elem_ne in mat.elems_ne:
			elem_ne.eps_ne_k = elem_ne.eps_ne_old.clone()
			elem_ne.eps_ne_rate = elem_ne.eps_ne_rate_old.clone()

		if len(mat.elems_ne) > 0:
			# Non-linear loop
			tol = 1e-10
			error = 2*tol
			ite = 0

			while error > tol:
				eps_k = []
				for elem_ne in mat.elems_ne:
					eps_k.append(elem_ne.eps_ne_k.flatten())
				eps_k = to.hstack(eps_k)

				eps = []
				for elem_ne in mat.elems_ne:
					elem_ne.compute_eps_ne_rate(stress, phi1, T)
					elem_ne.compute_G_B(stress, dt, theta, T)
					elem_ne.eps_ne_rate -= elem_ne.B
					elem_ne.increment_internal_variables(stress, stress, dt)

					eps.append(elem_ne.eps_ne_k.flatten())
				eps = to.hstack(eps)

				error = to.linalg.norm(eps - eps_k)
				ite += 1

		# Calculate eps_ne_rate
		eps_rate = to.zeros((n_elems, 3, 3))
		for elem_ne in mat.elems_ne:
			eps_rate += elem_ne.eps_ne_rate

		# Initialize total strain tensor
		eps_tot = to.zeros((n_elems, 3, 3))

		# Add elastic strains to total strain tensor
		for elem_e in mat.elems_e:
			eps_tot += elem_e.eps_e

		# Compute non-elastic strains
		for elem_ne in mat.elems_ne:
			elem_ne.compute_eps_ne_k(phi1, phi2)

		# Add non-elastic strains to total strain tensor
		for elem_ne in mat.elems_ne:
			eps_tot += elem_ne.eps_ne_k
			elem_ne.update_internal_variables()
			elem_ne.update_eps_ne_old(stress, stress, phi2)
			elem_ne.update_eps_ne_rate_old()

		print(eps_tot)

		eps_tot_list[i] = eps_tot

	# print(eps_tot_list.shape)
	# print(eps_tot_list)

	e_xx = -100*eps_tot_list[:,0,0,0][:]
	e_yy = -100*eps_tot_list[:,0,1,1][:]
	e_zz = -100*eps_tot_list[:,0,2,2][:]
	e_eq = e_zz - e_xx

	time_h = time_list[:]/ut.hour

	# Plot loading schedule
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3))
	fig.subplots_adjust(top=0.970, bottom=0.155, left=0.062, right=0.980, hspace=0.35, wspace=0.312)
	fig.patch.set_alpha(0.0)

	ax1.plot(time_h, e_eq, ".-")
	ax1.set_xlabel("Time (hour)")
	ax1.set_ylabel("Equivalent strain (%)")
	ax1.grid(True)

	ax2.plot(time_h, e_zz, ".-", color="steelblue")
	ax2.plot(time_h, e_xx, ".-", color="lightcoral")
	ax2.set_xlabel("Time (hour)")
	ax2.set_ylabel("Total strain (%)")
	ax2.grid(True)

	plt.show()




if __name__ == '__main__':
	main()