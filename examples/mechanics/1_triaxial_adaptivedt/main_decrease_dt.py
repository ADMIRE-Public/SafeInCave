import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
import dolfinx as do
import os
import torch as to


class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.eps_ve = do.fem.Function(self.DG0_3x3)
        self.eps_cr = do.fem.Function(self.DG0_3x3)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        self.eps_ve.x.array[:] = to.flatten(self.mat.elems_ne[0].eps_ne_k)
        self.eps_cr.x.array[:] = to.flatten(self.mat.elems_ne[1].eps_ne_k)
        self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[2].eps_ne_k)
        self.Fvp.x.array[:] = self.mat.elems_ne[2].Fvp


class LinearMomentumMixedMod(sf.LinearMomentumMixed):
    def __init__(self, grid, theta, stab_scaling=1.0):
        super().__init__(grid, theta, stab_scaling)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.eps_ve = do.fem.Function(self.DG0_3x3)
        self.eps_cr = do.fem.Function(self.DG0_3x3)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        self.eps_ve.x.array[:] = to.flatten(self.mat.elems_ne[0].eps_ne_k)
        self.eps_cr.x.array[:] = to.flatten(self.mat.elems_ne[1].eps_ne_k)
        self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[2].eps_ne_k)
        self.Fvp.x.array[:] = self.mat.elems_ne[2].Fvp


def run(formulation):
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cube")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Define output folder
    output_folder = os.path.join("output", "case_0", f"{formulation}_decrease_dt")

    # Time settings: this example intentionally shrinks dt.
    unit = "hour"
    t_0 = 0.0
    dt = 6.0
    t_final = 24.0
    t_control = sf.TimeControllerAdaptive(
        initial_dt=dt,
        max_dt=12.0,
        initial_time=t_0,
        final_time=t_final,
        time_unit=unit,
        growth_factor=1.5,
        shrink_factor=0.5,
        easy_ratio_threshold=0.25,
        hard_ratio_threshold=0.50,
        max_bisections=10,
        maxiter=3,
    )

    # Define momentum equation
    theta = 0.5
    if formulation == "P1":
        mom_eq = LinearMomentumMod(grid, theta=theta)
    elif formulation == "P1P1":
        mom_eq = LinearMomentumMixedMod(grid, theta=theta, stab_scaling=0.0)
    elif formulation == "P1P1_Stab":
        mom_eq = LinearMomentumMixedMod(grid, theta=theta, stab_scaling=1.0)

    # Define material properties
    mat = sf.Material(grid.n_elems)

    # Set material density
    rho = 2000.0 * to.ones(grid.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Constitutive model tuned for stronger nonlinear response.
    E = 45e9 * to.ones(grid.n_elems)
    nu = 0.3 * to.ones(grid.n_elems)
    spring_0 = sf.Spring(E, nu, "spring")

    # Create Kelvin-Voigt viscoelastic element
    eta = 2.5e11 * to.ones(grid.n_elems)
    E = 4e9 * to.ones(grid.n_elems)
    nu = 0.32 * to.ones(grid.n_elems)
    kelvin = sf.Viscoelastic(eta, E, nu, "kelvin")

    # Create creep
    A = 8.0e-12 * to.ones(grid.n_elems)
    Q = 28000 * to.ones(grid.n_elems)
    n = 1.5 * to.ones(grid.n_elems)
    creep_0 = sf.DislocationCreep(A, Q, n, "creep")

    # Create Desai's viscoplastic model
    mu_1 = 2.0e-6 * to.ones(grid.n_elems)
    N_1 = 3.8 * to.ones(grid.n_elems)
    n = 1.5 * to.ones(grid.n_elems)
    a_1 = 6.0e-05 * to.ones(grid.n_elems)
    eta = 0.55 * to.ones(grid.n_elems)
    beta_1 = 0.0100 * to.ones(grid.n_elems)
    beta = 0.995 * to.ones(grid.n_elems)
    m = -0.5 * to.ones(grid.n_elems)
    gamma = 0.095 * to.ones(grid.n_elems)
    alpha_0 = 0.030 * to.ones(grid.n_elems)
    sigma_t = 0.3 * to.ones(grid.n_elems)
    desai = sf.ViscoplasticDesai(
        mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai"
    )

    # Create constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(kelvin)
    mat.add_to_non_elastic(creep_0)
    mat.add_to_non_elastic(desai)

    # Set constitutive model
    mom_eq.set_material(mat)

    # Set body forces
    g_vec = [0.0, 0.0, 0.0]
    mom_eq.build_body_force(g_vec)

    # Set initial temperature field
    T0_field = 293 * to.ones(grid.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Boundary conditions
    bc_west = momBC.DirichletBC(
        boundary_name="WEST",
        component=0,
        values=[0.0, 0.0],
        time_values=[0.0, t_control.t_final],
    )

    bc_bottom = momBC.DirichletBC(
        boundary_name="BOTTOM",
        component=2,
        values=[0.0, 0.0],
        time_values=[0.0, t_control.t_final],
    )

    bc_south = momBC.DirichletBC(
        boundary_name="SOUTH",
        component=1,
        values=[0.0, 0.0],
        time_values=[0.0, t_control.t_final],
    )

    bc_east = momBC.NeumannBC(
        boundary_name="EAST",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[4.0 * ut.MPa, 4.0 * ut.MPa],
        time_values=[0.0, t_control.t_final],
        g=g_vec[2],
    )

    bc_north = momBC.NeumannBC(
        boundary_name="NORTH",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[4.0 * ut.MPa, 4.0 * ut.MPa],
        time_values=[0.0, t_control.t_final],
        g=g_vec[2],
    )

    bc_top = momBC.NeumannBC(
        boundary_name="TOP",
        direction=2,
        density=0.0,
        ref_pos=0.0,
        values=[4.1 * ut.MPa, 16 * ut.MPa, 16 * ut.MPa, 6 * ut.MPa, 6 * ut.MPa],
        time_values=[
            0 * ut.hour,
            2 * ut.hour,
            4 * ut.hour,
            5 * ut.hour,
            6 * ut.hour,
        ],
        g=g_vec[2],
    )

    bc_handler = momBC.BcHandler()
    bc_handler.add_boundary_condition(bc_west)
    bc_handler.add_boundary_condition(bc_bottom)
    bc_handler.add_boundary_condition(bc_south)
    bc_handler.add_boundary_condition(bc_east)
    bc_handler.add_boundary_condition(bc_north)
    bc_handler.add_boundary_condition(bc_top)

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_handler)

    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("eps_ve", "Viscoelastic strain (-)")
    output_mom.add_output_field("eps_cr", "Creep strain (-)")
    output_mom.add_output_field("eps_vp", "Viscoplastic strain (-)")
    output_mom.add_output_field("Fvp", "Yield function (-)")
    outputs = [output_mom]

    convergence_criterion = "strain_based"
    sim = sf.Simulator_M(
        mom_eq,
        t_control,
        outputs,
        compute_elastic_response=True,
        convergence_criterion=convergence_criterion,
    )
    sim.run()


def main():
    run("P1")


if __name__ == "__main__":
    main()
