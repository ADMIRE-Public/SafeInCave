from safeincave.CavernBC import (
    CavernHandler,
    Cavern_PT,
    Cavern_T,
    Cavern_MassFlux,
    CavernVolumeComputer,
)
from safeincave.Utils import create_field_nodes
import CoolProp.CoolProp as CP
from safeincave import GridHandlerGMSH
import dolfinx as do
import numpy as np
import os
import unittest


class Test_CavernVolumeComputer(unittest.TestCase):
    def setUp(self):
        test_dir = os.path.dirname(__file__)
        self.grid = GridHandlerGMSH("geom", os.path.join(test_dir, "files", "cube_caverns"))

        self.__calculate_expected_volume()
        self.__expected_normals()

    def __calculate_expected_volume(self):
        Lz = 1.0
        Ly = 2.0
        Lx = 1.0
        height = Lz / 2
        R = 0.1 * (Lx + Ly) / 2
        self.volume_expected = 2 * R * 2 * R * height

    def __expected_normals(self):
        self.normals_full = np.array(
            [
                [0.00000000e00, -1.00000000e00, 0.00000000e00],
                [0.00000000e00, -1.00000000e00, 0.00000000e00],
                [0.00000000e00, 0.00000000e00, 1.00000000e00],
                [0.00000000e00, 0.00000000e00, 1.00000000e00],
                [-1.00000000e00, -7.40148699e-16, 0.00000000e00],
                [-1.00000000e00, 0.00000000e00, 4.44089212e-16],
                [-0.00000000e00, -1.00000000e00, 0.00000000e00],
                [0.00000000e00, -1.00000000e00, 0.00000000e00],
                [0.00000000e00, 0.00000000e00, 1.00000000e00],
                [1.00000000e00, 0.00000000e00, -0.00000000e00],
                [1.00000000e00, 0.00000000e00, 0.00000000e00],
                [0.00000000e00, 0.00000000e00, 1.00000000e00],
                [-0.00000000e00, 0.00000000e00, -1.00000000e00],
                [-1.00000000e00, -0.00000000e00, -4.44089208e-16],
                [0.00000000e00, -0.00000000e00, -1.00000000e00],
                [1.00000000e00, 0.00000000e00, 0.00000000e00],
                [0.00000000e00, 0.00000000e00, -1.00000000e00],
                [-1.00000000e00, 7.40148667e-16, 0.00000000e00],
                [0.00000000e00, 1.00000000e00, -4.44089213e-16],
                [7.40148683e-16, 1.00000000e00, -0.00000000e00],
                [1.00000000e00, 0.00000000e00, 0.00000000e00],
                [-7.40148683e-16, 1.00000000e00, 0.00000000e00],
                [0.00000000e00, -0.00000000e00, -1.00000000e00],
                [0.00000000e00, 1.00000000e00, 4.44089207e-16],
            ]
        )
        self.normals_half = np.array(
            [
                [-0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, -0.0],
                [1.0, 0.0, 0.0],
                [-0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [-1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 0.0, -1.0],
                [-0.0, 0.0, -1.0],
                [0.0, 0.0, -1.0],
                [0.0, -0.0, -1.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, -0.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        self.normals_quarter = np.array(
            [
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 1.0],
                [0.0, -0.0, 1.0],
                [-0.0, 0.0, 1.0],
                [0.0, -1.0, 0.0],
                [0.0, -1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [-0.0, 0.0, -1.0],
                [0.0, -0.0, -1.0],
                [0.0, 0.0, -1.0],
                [0.0, 0.0, -1.0],
            ]
        )

    def test_compute_volume(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        volume_full = cvc.compute()
        self.assertAlmostEqual(volume_full, self.volume_expected, delta=1e-8)

        cvc_half = CavernVolumeComputer(
            self.grid, "Cavern_half", [10.0, 2.0, -780.0], sym_scale=2
        )
        volume_half = cvc_half.compute()
        self.assertAlmostEqual(volume_half, self.volume_expected, delta=1e-8)

        cvc_quarter = CavernVolumeComputer(
            self.grid, "Cavern_quarter", [0.0, 0.0, 0.5], sym_scale=4
        )
        volume_quarter = cvc_quarter.compute()
        self.assertAlmostEqual(volume_quarter, self.volume_expected, delta=1e-8)

    def test_normals(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        normals = cvc.calculate_normals()
        np.testing.assert_allclose(normals, self.normals_full, rtol=1e-3)

        cvc = CavernVolumeComputer(self.grid, "Cavern_half")
        normals = cvc.calculate_normals()
        np.testing.assert_allclose(normals, self.normals_half, rtol=1e-3)

        cvc = CavernVolumeComputer(self.grid, "Cavern_quarter")
        normals = cvc.calculate_normals()
        np.testing.assert_allclose(normals, self.normals_quarter, rtol=1e-3)

    def test_compute_volume_with_displacement(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        CG1_3x1 = do.fem.functionspace(
            cvc.grid.mesh, ("Lagrange", 1, (self.grid.domain_dim,))
        )
        u = do.fem.Function(CG1_3x1)
        u_array = u.x.array.reshape(-1, 3).copy()

        def fun_ux(x, y, z):
            return 0.5 - x * 0.2

        u_array[:, 0] = create_field_nodes(cvc.grid, fun_ux)
        u.x.array[:] = u_array.flatten()
        volume_computed = cvc.compute(u)

        R = 0.15
        p1 = 0.5 - R
        p2 = 0.5 + R
        u1 = 0.5 - p1 * 0.2
        u2 = 0.5 - p2 * 0.2
        D = (p2 + u2) - (p1 + u1)
        height = 0.75 - 0.25
        V_expected = D * 2 * R * height

        self.assertAlmostEqual(volume_computed, V_expected, delta=1e-8)

    def test_midpoint(self):
        cvc = CavernVolumeComputer(self.grid, "Cavern_full")
        midpoint = cvc.calculate_cavern_midpoint(2)
        self.assertEqual(midpoint, 0.5)

        cvc = CavernVolumeComputer(self.grid, "Cavern_half")
        midpoint = cvc.calculate_cavern_midpoint(2)
        self.assertEqual(midpoint, 0.5)

        cvc = CavernVolumeComputer(self.grid, "Cavern_quarter")
        midpoint = cvc.calculate_cavern_midpoint(2)
        self.assertEqual(midpoint, 0.5)


class Test_T(unittest.TestCase):
    def setUp(self):
        self.cavern_name = "Cave-1"
        self.cavern_1 = Cavern_T(self.cavern_name, [303, 400.0], [0.0, 10.0])

    def test_initialize(self):
        self.assertEqual(self.cavern_1.T, 303)


class Test_PT(unittest.TestCase):
    def setUp(self):
        test_dir = os.path.dirname(__file__)
        self.grid = GridHandlerGMSH("geom", os.path.join(test_dir, "files", "cube_caverns"))

        MPa = 1e6

        self.fluid_1 = "Water"
        self.cave_name_1 = "Cavern_full"
        self.cave_1 = Cavern_PT(
            grid=self.grid,
            cavern_name=self.cave_name_1,
            sym_scale=1,
            reference_point=None,
            fluid=self.fluid_1,
            P_values=[0 * MPa, 0 * MPa],
            T_values=[303, 400.0],
            time_values=[0.0, 10.0],
        )

        self.fluid_2 = "Methane"
        self.cave_name_2 = "Cavern_quarter"
        self.cave_2 = Cavern_PT(
            grid=self.grid,
            cavern_name=self.cave_name_2,
            sym_scale=4,
            reference_point=[0.0, 0.0, 0.5],
            fluid=self.fluid_2,
            P_values=[1 * MPa, 1 * MPa],
            T_values=[303, 400.0],
            time_values=[0.0, 10.0],
        )

        self.cavern_set = CavernHandler()
        self.cavern_set.add_cavern(self.cave_1)
        self.cavern_set.add_cavern(self.cave_2)

    def test_initialize(self):
        self.assertEqual(self.cave_1.AS.phase(), 0)  # CP.iphase_liquid
        self.assertEqual(self.cave_1.cavern_name, self.cave_name_1)
        self.assertEqual(self.cave_1.fluid, self.fluid_1)
        self.assertEqual(self.cave_1.AS.phase(), CP.get_phase_index("phase_liquid"))
        self.assertEqual(self.cave_1.density, 995.6946644467557)

        self.assertEqual(self.cave_2.AS.phase(), 2)  # CP.iphase_liquid
        self.assertEqual(self.cave_2.cavern_name, self.cave_name_2)
        self.assertEqual(self.cave_2.fluid, self.fluid_2)
        self.assertEqual(
            self.cave_2.AS.phase(), CP.get_phase_index("phase_supercritical_gas")
        )
        self.assertEqual(self.cave_2.density, 7.140352655623106)

    def test_update_cavern(self):
        self.cave_1.update_cavern(t=5.0, dt=None)
        self.assertEqual(self.cave_1.AS.phase(), CP.get_phase_index("phase_liquid"))
        self.assertEqual(self.cave_1.density, 972.8118876018121)

        self.cave_1.update_cavern(t=10.0, dt=None)
        self.assertEqual(self.cave_1.AS.phase(), CP.get_phase_index("phase_gas"))
        self.assertEqual(self.cave_1.density, 0.5549439034904987)

    def test_cavern_handler(self):
        self.assertEqual(len(self.cavern_set.caverns_PT), 2)
        self.cavern_set.update_caverns(2.0)

        self.assertEqual(
            self.cavern_set.caverns_PT[0].AS.phase(), CP.get_phase_index("phase_liquid")
        )
        self.assertEqual(self.cavern_set.caverns_PT[0].density, 988.3724191374077)

        self.assertEqual(
            self.cavern_set.caverns_PT[1].AS.phase(),
            CP.get_phase_index("phase_supercritical_gas"),
        )
        self.assertEqual(self.cavern_set.caverns_PT[1].density, 6.684645241990176)


class Test_MFlux(unittest.TestCase):
    def setUp(self):
        test_dir = os.path.dirname(__file__)
        self.grid = GridHandlerGMSH("geom", os.path.join(test_dir, "files", "cube_caverns"))

        CG1_3x1 = do.fem.functionspace(
            self.grid.mesh, ("Lagrange", 1, (self.grid.domain_dim,))
        )
        self.u = do.fem.Function(CG1_3x1)
        u_array = self.u.x.array.reshape(-1, 3).copy()

        def fun_ux(x, y, z):
            return 0.0

        u_array[:, 0] = create_field_nodes(self.grid, fun_ux)
        self.u.x.array[:] = u_array.flatten()

        self.MPa = 1e6
        minute = 60
        self.hour = 60 * minute

        t0 = 0.0
        tf = 10.0

        self.fluid_1 = "Water"
        self.cave_name_1 = "Cavern_full"
        self.cave_1 = Cavern_MassFlux(
            grid=self.grid,
            cavern_name=self.cave_name_1,
            sym_scale=1,
            reference_point=None,
            fluid=self.fluid_1,
            P_init=1 * self.MPa,
            T_init=303,
            T_in=303,
            Mflux_values=[0, 0.01, 0.01, 0.0, 0.0],
            time_values=np.linspace(t0, tf, 5),
        )

        self.fluid_2 = "Methane"
        self.cave_name_2 = "Cavern_half"
        self.cave_2 = Cavern_MassFlux(
            grid=self.grid,
            cavern_name=self.cave_name_2,
            sym_scale=2,
            reference_point=[10.0, 2.0, -780.0],
            fluid=self.fluid_2,
            P_init=1 * self.MPa,
            T_init=303,
            T_in=293,
            Mflux_values=[0, 0.01, 0.01, 0.0, 0.0],
            time_values=np.linspace(t0, tf, 5),
        )

        self.fluid_3 = "Hydrogen"
        self.cave_name_3 = "Cavern_quarter"
        self.cave_3 = Cavern_MassFlux(
            grid=self.grid,
            cavern_name=self.cave_name_3,
            sym_scale=4,
            reference_point=[0.0, 0.0, 0.5],
            fluid=self.fluid_3,
            P_init=1 * self.MPa,
            T_init=303,
            T_in=293,
            Mflux_values=[0, 0.01, 0.01, 0.0, 0.0],
            time_values=np.linspace(t0, tf, 5),
        )

        self.cavern_set = CavernHandler()
        self.cavern_set.add_cavern(self.cave_1)
        self.cavern_set.add_cavern(self.cave_2)
        self.cavern_set.add_cavern(self.cave_3)

        self.expected_density_1 = np.array(
            [
                996.140076,
                996.162298,
                996.206742,
                996.273409,
                996.362298,
                996.473409,
                996.58452,
                996.695631,
                996.806742,
                996.917854,
                997.028965,
                997.117854,
                997.18452,
                997.228965,
                997.251187,
                997.251187,
                997.251187,
                997.251187,
                997.251187,
            ]
        )
        self.expected_M_1 = np.array(
            [
                44.826303,
                44.827303,
                44.829303,
                44.832303,
                44.836303,
                44.841303,
                44.846303,
                44.851303,
                44.856303,
                44.861303,
                44.866303,
                44.870303,
                44.873303,
                44.875303,
                44.876303,
                44.876303,
                44.876303,
                44.876303,
                44.876303,
            ]
        )
        self.expected_P_1 = np.array(
            [
                1.0,
                1.050706,
                1.152135,
                1.304323,
                1.507324,
                1.761208,
                2.015242,
                2.269424,
                2.523755,
                2.778235,
                3.032864,
                3.236675,
                3.389595,
                3.491572,
                3.54257,
                3.54257,
                3.54257,
                3.54257,
                3.54257,
            ]
        )
        self.expected_T_1 = np.array(
            [
                303.0,
                303.001118,
                303.003356,
                303.006715,
                303.011197,
                303.016808,
                303.022426,
                303.028052,
                303.033685,
                303.039327,
                303.044976,
                303.049501,
                303.052898,
                303.055164,
                303.056297,
                303.056297,
                303.056297,
                303.056297,
                303.056297,
            ]
        )

        self.expected_density_2 = np.array(
            [
                7.140353,
                7.162575,
                7.207019,
                7.273686,
                7.362575,
                7.473686,
                7.584797,
                7.695908,
                7.807019,
                7.91813,
                8.029242,
                8.11813,
                8.184797,
                8.229242,
                8.251464,
                8.251464,
                8.251464,
                8.251464,
                8.251464,
            ]
        )
        self.expected_M_2 = np.array(
            [
                0.321316,
                0.322316,
                0.324316,
                0.327316,
                0.331316,
                0.336316,
                0.341316,
                0.346316,
                0.351316,
                0.356316,
                0.361316,
                0.365316,
                0.368316,
                0.370316,
                0.371316,
                0.371316,
                0.371316,
                0.371316,
                0.371316,
            ]
        )
        self.expected_P_2 = np.array(
            [
                1.0,
                1.004321,
                1.012962,
                1.025919,
                1.043187,
                1.064762,
                1.086324,
                1.107874,
                1.129412,
                1.150939,
                1.172454,
                1.189658,
                1.202557,
                1.211154,
                1.215452,
                1.215452,
                1.215452,
                1.215452,
                1.215452,
            ]
        )
        self.expected_T_2 = np.array(
            [
                303.0,
                303.249293,
                303.743011,
                304.471648,
                305.421513,
                306.575347,
                307.69356,
                308.777767,
                309.829487,
                310.85015,
                311.841104,
                312.613334,
                313.180937,
                313.55397,
                313.738903,
                313.738903,
                313.738903,
                313.738903,
                313.738903,
            ]
        )

        self.expected_density_3 = np.array(
            [
                0.875644,
                0.897866,
                0.942311,
                1.008977,
                1.097866,
                1.208977,
                1.320088,
                1.431199,
                1.542311,
                1.653422,
                1.764533,
                1.853422,
                1.920088,
                1.964533,
                1.986755,
                1.986755,
                1.986755,
                1.986755,
                1.986755,
            ]
        )
        self.expected_M_3 = np.array(
            [
                0.039404,
                0.040404,
                0.042404,
                0.045404,
                0.049404,
                0.054404,
                0.059404,
                0.064404,
                0.069404,
                0.074404,
                0.079404,
                0.083404,
                0.086404,
                0.088404,
                0.089404,
                0.089404,
                0.089404,
                0.089404,
                0.089404,
            ]
        )
        self.expected_P_3 = np.array(
            [
                1.0,
                1.038259,
                1.114802,
                1.229689,
                1.383032,
                1.575001,
                1.767319,
                1.960005,
                2.15307,
                2.346523,
                2.540369,
                2.695732,
                2.812422,
                2.890296,
                2.929258,
                2.929258,
                2.929258,
                2.929258,
                2.929258,
            ]
        )
        self.expected_T_3 = np.array(
            [
                303.0,
                305.7095,
                310.740987,
                317.450249,
                325.121802,
                333.121655,
                339.776256,
                345.401757,
                350.222219,
                354.400945,
                358.059806,
                360.675718,
                362.481451,
                363.618494,
                364.168332,
                364.168332,
                364.168332,
                364.168332,
                364.168332,
            ]
        )

    def test_results(self):
        t = 0.0
        dt = 0.5
        for _ in range(1, 20):
            # Advance time
            t += dt

            # Update cavern states at time t
            self.cavern_set.calculate_volumes(self.u)
            self.cavern_set.record_cavern_data(t)
            self.cavern_set.calculate_initial_conditions()
            self.cavern_set.update_caverns(t=t, dt=dt)

        # Assert history values are close to expected (pre-computed with CoolProp)
        np.testing.assert_allclose(
            self.cave_1.density_hist, self.expected_density_1, rtol=1e-5
        )
        np.testing.assert_allclose(self.cave_1.M_hist, self.expected_M_1, rtol=1e-5)
        np.testing.assert_allclose(
            self.cave_1.P_hist, self.expected_P_1 * self.MPa, rtol=1e-5
        )
        np.testing.assert_allclose(self.cave_1.T_hist, self.expected_T_1, rtol=1e-5)

        np.testing.assert_allclose(
            self.cave_2.density_hist, self.expected_density_2, rtol=1e-5
        )
        np.testing.assert_allclose(self.cave_2.M_hist, self.expected_M_2, rtol=1e-5)
        np.testing.assert_allclose(
            self.cave_2.P_hist, self.expected_P_2 * self.MPa, rtol=1e-5
        )
        np.testing.assert_allclose(self.cave_2.T_hist, self.expected_T_2, rtol=1e-5)

        np.testing.assert_allclose(
            self.cave_3.density_hist, self.expected_density_3, rtol=1e-5
        )
        np.testing.assert_allclose(self.cave_3.M_hist, self.expected_M_3, rtol=1e-5)
        np.testing.assert_allclose(
            self.cave_3.P_hist, self.expected_P_3 * self.MPa, rtol=1e-5
        )
        np.testing.assert_allclose(self.cave_3.T_hist, self.expected_T_3, rtol=1e-5)

        # import pandas as pd

        # # Create a DataFrame from the recorded data
        # df_1 = pd.DataFrame({
        #     "Time (s)": self.cave_1.t_hist,
        #     "Volume (m^3)": self.cave_1.V_hist,
        #     "Density (kg/m^3)": self.cave_1.density_hist,
        #     "Mass (kg)": self.cave_1.M_hist,
        #     "Pressure (MPa)": [P / self.MPa for P in self.cave_1.P_hist],
        #     "Temperature (K)": self.cave_1.T_hist
        # })
        # print(df_1)
        # print()

        # df_2 = pd.DataFrame({
        #     "Time (s)": self.cave_2.t_hist,
        #     "Volume (m^3)": self.cave_2.V_hist,
        #     "Density (kg/m^3)": self.cave_2.density_hist,
        #     "Mass (kg)": self.cave_2.M_hist,
        #     "Pressure (MPa)": [P / self.MPa for P in self.cave_2.P_hist],
        #     "Temperature (K)": self.cave_2.T_hist
        # })
        # print(df_2)
        # print()

        # df_3 = pd.DataFrame({
        #     "Time (s)": self.cave_3.t_hist,
        #     "Volume (m^3)": self.cave_3.V_hist,
        #     "Density (kg/m^3)": self.cave_3.density_hist,
        #     "Mass (kg)": self.cave_3.M_hist,
        #     "Pressure (MPa)": [P / self.MPa for P in self.cave_3.P_hist],
        #     "Temperature (K)": self.cave_3.T_hist
        # })
        # print(df_3)
