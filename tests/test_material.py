# import unittest
# import os
# import sys
# sys.path.append(os.path.join("..", "safeincave"))
# import torch as to
# import numpy as np
# from MaterialProps import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai

import safeincave as sc
import torch as to
import numpy as np
import unittest

class TestSpring(unittest.TestCase):
	def setUp(self):
		self.n_elems = 2
		E = to.tensor(102e9*np.ones(self.n_elems))
		nu = to.tensor(0.3*np.ones(self.n_elems))
		self.elem = sc.Spring(E, nu, name="spring")
		self.elem.initialize()

		self.stress = 1e6*to.tensor([[  [1., 4., 5.],
	                                 	[4., 2., 6.],
	                                 	[5., 6., 3.] ],
	                                  [ [6., 1., 2.],
	                                 	[1., 5., 3.],
	                                 	[2., 3., 4.] ]], dtype=to.float64)

		self.true_eps_e = to.tensor([[	 [-4.9020e-06,  5.0980e-05,  6.3725e-05],
								         [ 5.0980e-05,  7.8431e-06,  7.6471e-05],
								         [ 6.3725e-05,  7.6471e-05,  2.0588e-05]],

								        [[ 3.2353e-05,  1.2745e-05,  2.5490e-05],
								         [ 1.2745e-05,  1.9608e-05,  3.8235e-05],
								         [ 2.5490e-05,  3.8235e-05,  6.8627e-06]]], dtype=to.float64)

	def test_eps_e(self):
		self.elem.compute_eps_e(self.stress)
		to.testing.assert_close(self.elem.eps_e, self.true_eps_e, rtol=1e-6, atol=1e-9)

class TestViscoelastic(unittest.TestCase):
	def setUp(self):
		self.n_elems = 1
		props = {
			"E":   to.tensor(10e9*np.ones(self.n_elems)),
			"nu":  to.tensor(0.32*np.ones(self.n_elems)),
			"eta": to.tensor(105e11*np.ones(self.n_elems))
		}
		E = to.tensor(10e9*np.ones(self.n_elems))
		nu = to.tensor(0.32*np.ones(self.n_elems))
		eta = to.tensor(105e11*np.ones(self.n_elems))
		self.elem = sc.Viscoelastic(eta, E, nu)

		self.stress = 1e6*to.tensor([[  [1., 4., 5.],
	                                 	[4., 2., 6.],
	                                 	[5., 6., 3.]] ], dtype=to.float64)

		self.zeros3x3 = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.Temp = 298*to.ones(self.n_elems, dtype=to.float64)

		self.theta = 0.5
		self.dt = 7200.
		self.phi2 = (1 - self.theta)*self.dt
		self.phi1 = self.theta*self.dt

		self.true_G = to.tensor([[	 [ 2.0666e-14, -5.8081e-15, -5.8081e-15,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [-5.8081e-15,  2.0666e-14, -5.8081e-15,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [-5.8081e-15, -5.8081e-15,  2.0666e-14,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  2.6474e-14,  0.0000e+00, 0.0000e+00],
							         [-0.0000e+00, -0.0000e+00, -0.0000e+00, -0.0000e+00,  2.6474e-14, -0.0000e+00],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, 2.6474e-14]]], dtype=to.float64)

		self.true_eps_ve_rate = to.tensor([[ [-8.3746e-09,  1.0590e-07,  1.3237e-07],
									         [ 1.0590e-07,  1.8100e-08,  1.5884e-07],
									         [ 1.3237e-07,  1.5884e-07,  4.4574e-08]]], dtype=to.float64)

		self.true_eps_k = to.tensor([[ [-3.0148e-05,  3.8123e-04,  4.7653e-04],
								         [ 3.8123e-04,  6.5158e-05,  5.7184e-04],
								         [ 4.7653e-04,  5.7184e-04,  1.6047e-04]]], dtype=to.float64)

		self.true_eps_ve = to.tensor([[	 [-6.0297e-05,  7.6245e-04,  9.5307e-04],
								         [ 7.6245e-04,  1.3032e-04,  1.1437e-03],
								         [ 9.5307e-04,  1.1437e-03,  3.2093e-04]]], dtype=to.float64)


	def test_full(self):
		self.elem.compute_G_B(self.stress, self.dt, self.theta, self.Temp)
		to.testing.assert_close(self.elem.G, self.true_G, rtol=1e-18, atol=1e-18)

		self.elem.compute_eps_ne_rate(self.stress, self.phi1, self.Temp)
		to.testing.assert_close(self.elem.eps_ne_rate, self.true_eps_ve_rate, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_ne_k(self.phi1, self.phi2)
		to.testing.assert_close(self.elem.eps_ne_k, self.true_eps_k, rtol=1e-8, atol=1e-8)

		self.elem.update_eps_ne_old(self.stress, self.zeros3x3, self.phi2)
		to.testing.assert_close(self.elem.eps_ne_old, self.true_eps_ve, rtol=1e-7, atol=1e-7)

		self.elem.update_eps_ne_rate_old()
		to.testing.assert_close(self.elem.eps_ne_rate_old, self.true_eps_ve_rate, rtol=1e-7, atol=1e-7)


class TestDislocationCreep(unittest.TestCase):
	def setUp(self):
		self.n_elems = 1
		A = 1.9e-20*to.ones(self.n_elems)
		n = 3.0*to.ones(self.n_elems)
		Q = 51600*to.ones(self.n_elems)
		self.elem = sc.DislocationCreep(A, Q, n, name="creep")

		self.stress = 1e6*to.tensor([[  [1., 4., 5.],
	                                 	[4., 2., 6.],
	                                 	[5., 6., 3.]] ], dtype=to.float64)

		self.zeros3x3 = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.Temp = 298*to.ones(self.n_elems, dtype=to.float64)

		self.theta = 0.5
		self.dt = 7200.
		self.phi2 = (1 - self.theta)*self.dt
		self.phi1 = self.theta*self.dt

		self.true_G = to.tensor([[	[ 2.7650e-15, -1.3564e-15, -1.4086e-15, -8.3471e-16, -1.0434e-15, -1.2521e-15],
         							[-1.3564e-15,  2.7128e-15, -1.3564e-15,  0.0000e+00,  0.0000e+00, 0.0000e+00],
         							[-1.4086e-15, -1.3564e-15,  2.7650e-15,  8.3471e-16,  1.0434e-15, 1.2521e-15],
         							[-2.0868e-16,  0.0000e+00,  2.0868e-16,  1.1477e-14,  4.1735e-15, 5.0083e-15],
         							[-2.6085e-16,  0.0000e+00,  2.6085e-16,  4.1735e-15,  1.3355e-14, 6.2603e-15],
         							[-3.1302e-16,  0.0000e+00,  3.1302e-16,  5.0083e-15,  6.2603e-15, 1.5651e-14]]], dtype=to.float64)

		self.true_eps_ne_rate = to.tensor([[ [-4.0692e-09,  1.6277e-08,  2.0346e-08],
									         [ 1.6277e-08,  0.0000e+00,  2.4415e-08],
									         [ 2.0346e-08,  2.4415e-08,  4.0692e-09]]], dtype=to.float64)

		self.true_eps_ne_k = to.tensor([[ [-1.4649e-05,  5.8597e-05,  7.3246e-05],
								          [ 5.8597e-05,  0.0000e+00,  8.7895e-05],
								          [ 7.3246e-05,  8.7895e-05,  1.4649e-05]]], dtype=to.float64)

		self.true_eps_ne = to.tensor([[	 [-8.7519e-05,  4.0867e-04,  5.1084e-04],
								         [ 4.0867e-04,  1.3643e-12,  6.1301e-04],
								         [ 5.1084e-04,  6.1301e-04,  8.7519e-05]]], dtype=to.float64)

	def test_full(self):
		self.elem.compute_G_B(self.stress, self.dt, self.theta, self.Temp)
		to.testing.assert_close(self.elem.G, self.true_G, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_ne_rate(self.stress, self.phi1, self.Temp)
		to.testing.assert_close(self.elem.eps_ne_rate, self.true_eps_ne_rate, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_ne_k(self.phi1, self.phi2)
		to.testing.assert_close(self.elem.eps_ne_k, self.true_eps_ne_k, rtol=1e-8, atol=1e-8)

		self.elem.update_eps_ne_old(self.stress, self.zeros3x3, self.phi2)
		to.testing.assert_close(self.elem.eps_ne_old, self.true_eps_ne, rtol=1e-4, atol=1e-4)

		to.testing.assert_close(self.elem.eps_ne_rate_old, self.zeros3x3, rtol=1e-10, atol=1e-10)
		self.elem.update_eps_ne_rate_old()
		to.testing.assert_close(self.elem.eps_ne_rate_old, self.true_eps_ne_rate, rtol=1e-10, atol=1e-10)


class TestViscoplasticDesai(unittest.TestCase):
	def setUp(self):
		self.n_elems = 1
		mu_1 = 5.3665857009859815e-11*to.ones(self.n_elems, dtype=to.float64)
		N_1 = 3.1*to.ones(self.n_elems, dtype=to.float64)
		n = 3.0*to.ones(self.n_elems, dtype=to.float64)
		a_1 = 1.965018496922832e-05*to.ones(self.n_elems, dtype=to.float64)
		eta = 0.8275682807874163*to.ones(self.n_elems, dtype=to.float64)
		beta_1 = 0.0048*to.ones(self.n_elems, dtype=to.float64)
		beta = 0.995*to.ones(self.n_elems, dtype=to.float64)
		m = -0.5*to.ones(self.n_elems, dtype=to.float64)
		gamma = 0.095*to.ones(self.n_elems, dtype=to.float64)
		alpha_0 = 0.0022*to.ones(self.n_elems, dtype=to.float64)
		sigma_t = 5.0*to.ones(self.n_elems, dtype=to.float64)
		self.elem = sc.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, name="desai")

		self.stress = -1e7*to.tensor([[ [1., 0., 0.],
	                                 	[0., 1., 0.],
	                                 	[0., 0., 3.]] ], dtype=to.float64)

		self.zeros3x3 = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		self.Temp = 298*to.ones(self.n_elems, dtype=to.float64)

		self.theta = 0.5
		self.dt = 7200.
		self.phi2 = (1 - self.theta)*self.dt
		self.phi1 = self.theta*self.dt

		self.true_G = to.tensor([[	[-1.1784e-02, -1.1784e-02, -1.1784e-02, -2.3568e-02, -2.3568e-02, -2.3568e-02],
							        [-1.1784e-02, -1.1784e-02, -1.1784e-02, -2.3568e-02, -2.3568e-02, -2.3568e-02],
							        [ 1.0323e-01,  1.0323e-01,  1.0323e-01,  2.0645e-01,  2.0645e-01,  2.0645e-01],
							        [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  2.7239e-09,  0.0000e+00,  0.0000e+00],
							        [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  2.2995e-09,  0.0000e+00],
							        [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  2.2995e-09]]], dtype=to.float64)

		self.true_eps_ne_rate = to.tensor([[  [ 0.00117867, -0.,         -0.        ],
											  [-0.,          0.00117867, -0.        ],
											  [-0.,         -0.,         -0.01031885]]], dtype=to.float64)

		self.true_eps_ne_k = to.tensor([[   [  4.24319634,   0.,           0.        ],
								            [  0.,           4.24319634,   0.        ],
								            [  0.,           0.,         -37.14785876]]], dtype=to.float64)

		self.true_eps_ne = to.tensor([[   [ 2.12113385e+09,  0.00000000e+00,  0.00000000e+00],
										  [ 0.00000000e+00,  2.12113385e+09,  0.00000000e+00],
										  [ 0.00000000e+00,  0.00000000e+00, -1.85806882e+10]]], dtype=to.float64)

		self.true_Fvp = to.tensor([185.2260], dtype=to.float64)
		self.true_alpha = to.tensor([0.0022], dtype=to.float64)
		self.true_qsi = to.tensor([75.2588], dtype=to.float64)

	def test_full(self):
		self.elem.compute_G_B(self.stress, self.dt, self.theta, self.Temp)
		to.testing.assert_close(self.elem.G, self.true_G, rtol=1e-4, atol=1e-10)

		self.elem.compute_eps_ne_rate(self.stress, self.phi1, self.Temp)
		to.testing.assert_close(self.elem.eps_ne_rate, self.true_eps_ne_rate, rtol=1e-4, atol=1e-10)

		self.elem.compute_eps_ne_k(self.phi1, self.phi2)
		to.testing.assert_close(self.elem.eps_ne_k, self.true_eps_ne_k, rtol=1e-4, atol=1e-8)

		self.elem.update_eps_ne_old(self.stress, self.zeros3x3, self.phi2)
		to.testing.assert_close(self.elem.eps_ne_old, self.true_eps_ne, rtol=1e-4, atol=1e-4)

		to.testing.assert_close(self.elem.eps_ne_rate_old, self.zeros3x3, rtol=1e-10, atol=1e-10)

		to.testing.assert_close(self.elem.Fvp, self.true_Fvp, rtol=1e-3, atol=1e-4)
		to.testing.assert_close(self.elem.alpha, self.true_alpha, rtol=1e-3, atol=1e-4)
		to.testing.assert_close(self.elem.qsi, self.true_qsi, rtol=1e-3, atol=1e-4)

