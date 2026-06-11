from safeincave import CavernThermodynamics
import unittest


class Test_T(unittest.TestCase):
    def setUp(self):
        self.fluid = "Air"
        self.thermo = CavernThermodynamics(self.fluid)

    def test_rho_u_h(self):
        rho, u, h = self.thermo.rho_u_h(101325, 300)
        self.assertEqual(rho, 1.1769955883877592)
        self.assertEqual(u, 340209.941098455)
        self.assertEqual(h, 426297.7743916913)

    def test_solve_withdrawal(self):
        P0 = 101325
        T0 = 300
        V0 = 1000.0
        V1 = 990.0
        dm = -1.0
        Q = 10000.0
        P1, T1, rho1 = self.thermo.solve_withdrawal(
            P0=P0, T0=T0, V0=V0, m_prod=dm, Q=Q, V1=V1
        )
        self.assertAlmostEqual(P1, 102643.79952048359, places=6)
        self.assertAlmostEqual(T1, 301.1196554705438, places=6)
        self.assertAlmostEqual(rho1, 1.1878743317048073, places=6)

    def test_solve_injection(self):
        P0 = 101325
        T0 = 300
        V0 = 1000.0
        V1 = 1000.0
        T_in = 420
        dm = 10000.0
        Q = 100000000.0
        P1, T1, rho1 = self.thermo.solve_injection(
            P0=P0, T0=T0, V0=V0, Tin=T_in, m_inj=dm, Q=Q, V1=V1
        )
        self.assertAlmostEqual(P1, 1828353.1178692507, places=6)
        self.assertAlmostEqual(T1, 566.0800589800913, places=6)
        self.assertAlmostEqual(rho1, 11.17699558838776, places=6)
