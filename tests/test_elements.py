import unittest
import os
import sys
sys.path.append(os.path.join("..", "safeincave"))
import torch as to
import numpy as np
from Elements import Spring, Viscoelastic, DislocationCreep, ViscoplasticDesai

class TestSpring(unittest.TestCase):
	def setUp(self):
		self.n_elems = 2
		props = {
			"E":  to.tensor(102e9*np.ones(self.n_elems)),
			"nu": to.tensor(0.3*np.ones(self.n_elems))
		}
		self.elem = Spring(props)
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
		self.elem = Viscoelastic(props)

		self.stress = 1e6*to.tensor([[  [1., 4., 5.],
	                                 	[4., 2., 6.],
	                                 	[5., 6., 3.]] ], dtype=to.float64)

		self.zeros3x3 = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		theta = 0.5
		self.dt = 7200.
		self.phi2 = (1 - theta)*self.dt
		self.phi1 = theta*self.dt

		self.true_G = to.tensor([[	 [ 2.0666e-14, -5.8081e-15, -5.8081e-15,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [-5.8081e-15,  2.0666e-14, -5.8081e-15,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [-5.8081e-15, -5.8081e-15,  2.0666e-14,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  2.6474e-14,  0.0000e+00, 0.0000e+00],
							         [-0.0000e+00, -0.0000e+00, -0.0000e+00, -0.0000e+00,  2.6474e-14, -0.0000e+00],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, 2.6474e-14]]], dtype=to.float64)

		self.true_eps_ve_rate = to.tensor([[ [-8.3746e-09,  1.0590e-07,  1.3237e-07],
									         [ 1.0590e-07,  1.8100e-08,  1.5884e-07],
									         [ 1.3237e-07,  1.5884e-07,  4.4574e-08]]], dtype=to.float64)

		self.true_eps_bar = to.tensor([[ [-3.0148e-05,  3.8123e-04,  4.7653e-04],
								         [ 3.8123e-04,  6.5158e-05,  5.7184e-04],
								         [ 4.7653e-04,  5.7184e-04,  1.6047e-04]]], dtype=to.float64)

		self.true_eps_ve = to.tensor([[	 [-6.0297e-05,  7.6245e-04,  9.5307e-04],
								         [ 7.6245e-04,  1.3032e-04,  1.1437e-03],
								         [ 9.5307e-04,  1.1437e-03,  3.2093e-04]]], dtype=to.float64)


	def test_full(self):
		self.elem.compute_G_B(self.stress, self.phi2)
		to.testing.assert_close(self.elem.G, self.true_G, rtol=1e-18, atol=1e-18)

		self.elem.compute_eps_ve_rate(self.stress, self.phi1)
		to.testing.assert_close(self.elem.eps_ve_rate, self.true_eps_ve_rate, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_bar(self.phi1, self.phi2)
		to.testing.assert_close(self.elem.eps_bar, self.true_eps_bar, rtol=1e-8, atol=1e-8)

		self.elem.compute_eps_ve(self.stress, self.zeros3x3, self.phi2)
		to.testing.assert_close(self.elem.eps_ve, self.true_eps_ve, rtol=1e-7, atol=1e-7)

		to.testing.assert_close(self.elem.eps_ve_old, self.zeros3x3, rtol=1e-7, atol=1e-7)
		self.elem.update_eps_ve_old()
		to.testing.assert_close(self.elem.eps_ve_old, self.true_eps_ve, rtol=1e-7, atol=1e-7)

		to.testing.assert_close(self.elem.eps_ve_rate_old, self.zeros3x3, rtol=1e-10, atol=1e-10)
		self.elem.update_eps_ve_rate_old()
		to.testing.assert_close(self.elem.eps_ve_rate_old, self.true_eps_ve_rate, rtol=1e-10, atol=1e-10)


class TestDislocationCreep(unittest.TestCase):
	def setUp(self):
		self.n_elems = 1
		props = {
			"A": to.tensor(1.9e-20*np.ones(self.n_elems)),
			"n": to.tensor(3.0*np.ones(self.n_elems)),
			"T": to.tensor(298*np.ones(self.n_elems)),
			"Q": to.tensor(51600*np.ones(self.n_elems)),
			"R": to.tensor(8.32*np.ones(self.n_elems))
		}
		self.elem = DislocationCreep(props)

		self.stress = 1e6*to.tensor([[  [1., 4., 5.],
	                                 	[4., 2., 6.],
	                                 	[5., 6., 3.]] ], dtype=to.float64)

		self.zeros3x3 = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		theta = 0.5
		self.dt = 7200.
		self.phi2 = (1 - theta)*self.dt
		self.phi1 = theta*self.dt

		self.true_G = to.tensor([[	 [-4.0692e-07, -4.0692e-07, -4.0692e-07, -8.1384e-07, -8.1384e-07, -8.1384e-07],
							         [-1.3564e-15,  2.7128e-15, -1.3564e-15,  0.0000e+00,  0.0000e+00, 0.0000e+00],
							         [ 4.0692e-07,  4.0692e-07,  4.0692e-07,  8.1384e-07,  8.1384e-07, 8.1384e-07],
							         [ 1.6277e-06,  1.6277e-06,  1.6277e-06,  3.2554e-06,  3.2554e-06, 3.2554e-06],
							         [ 2.0346e-06,  2.0346e-06,  2.0346e-06,  4.0692e-06,  4.0692e-06, 4.0692e-06],
							         [ 2.4415e-06,  2.4415e-06,  2.4415e-06,  4.8830e-06,  4.8830e-06, 4.8830e-06]]], dtype=to.float64)

		self.true_eps_ie_rate = to.tensor([[ [-4.0692e-09,  1.6277e-08,  2.0346e-08],
									         [ 1.6277e-08,  0.0000e+00,  2.4415e-08],
									         [ 2.0346e-08,  2.4415e-08,  4.0692e-09]]], dtype=to.float64)

		self.true_eps_bar = to.tensor([[ [-1.4649e-05,  5.8597e-05,  7.3246e-05],
								         [ 5.8597e-05,  0.0000e+00,  8.7895e-05],
								         [ 7.3246e-05,  8.7895e-05,  1.4649e-05]]], dtype=to.float64)

		self.true_eps_ie = to.tensor([[	 [-5.2737e+04,  2.1095e+05,  2.6368e+05],
								         [ 2.1095e+05,  1.3631e-12,  3.1642e+05],
								         [ 2.6368e+05,  3.1642e+05,  5.2737e+04]]], dtype=to.float64)

	def test_full(self):
		self.elem.compute_G_B(self.stress, self.phi2)
		to.testing.assert_close(self.elem.G, self.true_G, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_ie_rate(self.stress)
		to.testing.assert_close(self.elem.eps_ie_rate, self.true_eps_ie_rate, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_bar(self.phi1, self.phi2)
		to.testing.assert_close(self.elem.eps_bar, self.true_eps_bar, rtol=1e-8, atol=1e-8)

		self.elem.compute_eps_ie(self.stress, self.zeros3x3, self.phi2)
		to.testing.assert_close(self.elem.eps_ie, self.true_eps_ie, rtol=1e-4, atol=1e-4)

		to.testing.assert_close(self.elem.eps_ie_old, self.zeros3x3, rtol=1e-4, atol=1e-4)
		self.elem.update_eps_ie_old()
		to.testing.assert_close(self.elem.eps_ie_old, self.true_eps_ie, rtol=1e-4, atol=1e-4)

		to.testing.assert_close(self.elem.eps_ie_rate_old, self.zeros3x3, rtol=1e-10, atol=1e-10)
		self.elem.update_eps_ie_rate_old()
		to.testing.assert_close(self.elem.eps_ie_rate_old, self.true_eps_ie_rate, rtol=1e-10, atol=1e-10)


class TestViscoplasticDesai(unittest.TestCase):
	def setUp(self):
		self.n_elems = 1
		props = {
			"mu_1": 	to.tensor(5.3665857009859815e-11*np.ones(self.n_elems)),
			"N_1": 		to.tensor(3.1*np.ones(self.n_elems)),
			"n": 		to.tensor(3.0*np.ones(self.n_elems)),
			"a_1":		to.tensor(1.965018496922832e-05*np.ones(self.n_elems)),
			"eta": 		to.tensor(0.8275682807874163*np.ones(self.n_elems)),
			"beta_1": 	to.tensor(0.0048*np.ones(self.n_elems)),
			"beta": 	to.tensor(0.995*np.ones(self.n_elems)),
			"m": 		to.tensor(-0.5*np.ones(self.n_elems)),
			"gamma": 	to.tensor(0.095*np.ones(self.n_elems)),
			"alpha_0": 	to.tensor(0.0022*np.ones(self.n_elems)),
			"k_v": 		to.tensor(0.0*np.ones(self.n_elems)),
			"sigma_t": 	to.tensor(5.0*np.ones(self.n_elems))
		}
		self.elem = ViscoplasticDesai(props)

		self.stress = 1e6*to.tensor([[  [1., 0., 0.],
	                                 	[0., 2., 0.],
	                                 	[0., 0., 3.]] ], dtype=to.float64)

		self.zeros3x3 = to.zeros((self.n_elems, 3, 3), dtype=to.float64)

		theta = 0.5
		self.dt = 7200.
		self.phi2 = (1 - theta)*self.dt
		self.phi1 = theta*self.dt

		self.true_G = to.tensor([[	 [-5.8281e-10, -5.8281e-10, -5.8281e-10, -1.1656e-09, -1.1656e-09, -1.1656e-09],
							         [ 2.9391e-10,  2.9391e-10,  2.9391e-10,  5.8782e-10,  5.8782e-10, 5.8782e-10],
							         [-3.0257e-10, -3.0257e-10, -3.0257e-10, -6.0513e-10, -6.0513e-10, -6.0513e-10],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  9.7748e-16,  0.0000e+00, 0.0000e+00],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  1.5622e-16, 0.0000e+00],
							         [ 0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, -6.6504e-16]]], dtype=to.float64)

		self.true_eps_ie_rate = to.tensor([[ [-1.6245e-10, -0.0000e+00, -0.0000e+00],
									         [-0.0000e+00,  8.1923e-11, -0.0000e+00],
									         [-0.0000e+00, -0.0000e+00, -8.4336e-11]]], dtype=to.float64)

		self.true_eps_bar = to.tensor([[ [-5.8481e-07,  0.0000e+00,  0.0000e+00],
								         [ 0.0000e+00,  2.9492e-07,  0.0000e+00],
								         [ 0.0000e+00,  0.0000e+00, -3.0361e-07]]], dtype=to.float64)

		self.true_eps_ie = to.tensor([[  [-12.5886,   0.0000,   0.0000],
								         [  0.0000,   6.3485,   0.0000],
								         [  0.0000,   0.0000,  -6.5354]]], dtype=to.float64)

		self.true_Fvp = to.tensor([0.9026], dtype=to.float64)
		self.true_alpha = to.tensor([0.0022], dtype=to.float64)
		self.true_qsi = to.tensor([7.2192e-07], dtype=to.float64)

	def test_full(self):
		self.elem.compute_G_B(self.stress, self.phi2)
		to.testing.assert_close(self.elem.G, self.true_G, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_ie_rate(self.stress)
		to.testing.assert_close(self.elem.eps_ie_rate, self.true_eps_ie_rate, rtol=1e-10, atol=1e-10)

		self.elem.compute_eps_bar(self.phi1, self.phi2)
		to.testing.assert_close(self.elem.eps_bar, self.true_eps_bar, rtol=1e-8, atol=1e-8)

		self.elem.compute_eps_ie(self.stress, self.zeros3x3, self.phi2)
		to.testing.assert_close(self.elem.eps_ie, self.true_eps_ie, rtol=1e-4, atol=1e-4)

		to.testing.assert_close(self.elem.eps_ie_old, self.zeros3x3, rtol=1e-4, atol=1e-4)
		self.elem.update_eps_ie_old()
		to.testing.assert_close(self.elem.eps_ie_old, self.true_eps_ie, rtol=1e-4, atol=1e-4)

		to.testing.assert_close(self.elem.eps_ie_rate_old, self.zeros3x3, rtol=1e-10, atol=1e-10)
		self.elem.update_eps_ie_rate_old()
		to.testing.assert_close(self.elem.eps_ie_rate_old, self.true_eps_ie_rate, rtol=1e-10, atol=1e-10)

		to.testing.assert_close(self.elem.Fvp, self.true_Fvp, rtol=1e-4, atol=1e-4)
		to.testing.assert_close(self.elem.alpha, self.true_alpha, rtol=1e-4, atol=1e-4)
		to.testing.assert_close(self.elem.qsi, self.true_qsi, rtol=1e-11, atol=1e-11)

