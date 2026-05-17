"""
Minimal working example: salt cavern with the Munson-Dawson constitutive model.

Model stack:
    Spring (linear elastic) + MunsonDawsonCreep


Two parameter sets are provided below. Switch with the SCENARIO flag:
    "B" (default):  cyclic-triaxial calibration at T=21 C, sigma_3=12 MPa
                    (Van den Brekel, 2026, thesis dataset).
    "A":            CCC Zuidwending dataset.

Geometry: grids/cavern_irregular. Operation phase only (30 days, dt=2 h),
with two pressure cycles between 10 MPa and 7 MPa. The initial elastic
stress is obtained from a one-shot solve_elastic_response() and used to
initialise MD's zeta in steady state.

Run: python main.py  -> outputs in ./output/munson_dawson_example/
"""

import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import os


SCENARIO = "B"   # "A" or "B"

sec_per_year = 365.25 * 24.0 * 3600.0


def build_elastic_material(mom_eq):
    """Spring-only material. Returns (mat, salt_density, E0, nu0)."""

    mat = sf.Material(mom_eq.n_elems)

    salt_density = 2000.0
    rho = salt_density * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    E0  = 20.425 * ut.GPa * to.ones(mom_eq.n_elems)
    nu0 = 0.25 * to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring_0)
    return mat, salt_density, E0, nu0


def build_munson_dawson(mom_eq, scenario, E0, nu0):
    """Construct the MunsonDawsonCreep element (not yet attached)."""

    if scenario == "A":
        # Scenario A: CCC Zuidwending
        nmd   = 4.99
        A_md  = (18.31 * (1e-6)**nmd / sec_per_year) * to.ones(mom_eq.n_elems, dtype=to.float64)
        Q_md  = (6356.0 * 8.32) * to.ones(mom_eq.n_elems)
        n_md  = nmd * to.ones(mom_eq.n_elems)
        K0_md = 7.0e-7  * to.ones(mom_eq.n_elems)
        c_md  = 9.02e-3 * to.ones(mom_eq.n_elems)
        m_md  = 3.0     * to.ones(mom_eq.n_elems)
        aw_md = -13.2   * to.ones(mom_eq.n_elems)
        bw_md = -7.738  * to.ones(mom_eq.n_elems)
        d_md  = 0.58    * to.ones(mom_eq.n_elems)
    elif scenario == "B":
        # Scenario B: cyclic triaxial at T=21 C, sigma_3=12 MPa
        nmd   = 5.6897
        A_md  = (17.28 * (1e-6)**nmd / sec_per_year) * to.ones(mom_eq.n_elems, dtype=to.float64)
        Q_md  = (6495.0 * 8.32) * to.ones(mom_eq.n_elems)
        n_md  = nmd * to.ones(mom_eq.n_elems)
        K0_md = 2253.87 * to.ones(mom_eq.n_elems)
        c_md  = 9.02e-3 * to.ones(mom_eq.n_elems)
        m_md  = 2.466   * to.ones(mom_eq.n_elems)
        aw_md = 179.70  * to.ones(mom_eq.n_elems)
        bw_md = 60.00   * to.ones(mom_eq.n_elems)
        d_md  = 299.95  * to.ones(mom_eq.n_elems)
    else:
        raise ValueError(f"Unknown SCENARIO: {scenario!r}. Choose 'A' or 'B'.")

    mu_md = E0 / (2.0 * (1.0 + nu0))
    return sf.MunsonDawsonCreep(
        A=A_md, Q=Q_md, n=n_md,
        K0=K0_md, c=c_md, m=m_md,
        alpha_w=aw_md, beta_w=bw_md, delta=d_md,
        mu=mu_md, name="munson_dawson",
    )


def initialize_md_steady_state(md, mom_eq, T_kelvin):
    """Set zeta = eps_t_star at the current stress so F = 1 at t = 0
    (MD analogue of Desai.compute_initial_hardening)."""

    stress = ut.numpy2torch(
        mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
    s_xx, s_yy, s_zz = stress[:, 0, 0], stress[:, 1, 1], stress[:, 2, 2]
    s_xy, s_xz, s_yz = stress[:, 0, 1], stress[:, 0, 2], stress[:, 1, 2]
    sigma_vm = to.sqrt(
        0.5 * (
            (s_xx - s_yy) ** 2 + (s_xx - s_zz) ** 2 + (s_yy - s_zz) ** 2
            + 6.0 * (s_xy ** 2 + s_xz ** 2 + s_yz ** 2)
        )
    )
    sigma_safe = to.clamp(sigma_vm, min=1.0)
    mu_safe    = to.clamp(md.mu, min=1.0)
    Temp       = T_kelvin * to.ones(mom_eq.n_elems, dtype=to.float64)

    ratio      = to.clamp(sigma_safe / mu_safe, min=1e-30)
    eps_t_star = md.K0 * to.exp(md.c * Temp) * (ratio ** md.m)
    eps_t_star = to.clamp(eps_t_star, min=1e-50)

    md.zeta     = eps_t_star.clone()
    md.zeta_old = md.zeta.clone()


def build_bcs(mom_eq, t_final, salt_density, g_vec,
              side_burden, over_burden,
              cavern_values, cavern_times, gas_density):
    """Cavern BCs: symmetry (W/S/B), lithostatic Neumann (E/N/Top),
    cavern-pressure Neumann from (cavern_values, cavern_times)."""

    bc_west = momBC.DirichletBC(
        boundary_name="West", component=0,
        values=[0.0, 0.0], time_values=[0.0, t_final])

    bc_south = momBC.DirichletBC(
        boundary_name="South", component=1,
        values=[0.0, 0.0], time_values=[0.0, t_final])

    bc_bottom = momBC.DirichletBC(
        boundary_name="Bottom", component=2,
        values=[0.0, 0.0], time_values=[0.0, t_final])

    bc_east = momBC.NeumannBC(
        boundary_name="East", direction=2,
        density=salt_density, ref_pos=660.0,
        values=[side_burden, side_burden],
        time_values=[0.0, t_final], g=g_vec[2])

    bc_north = momBC.NeumannBC(
        boundary_name="North", direction=2,
        density=salt_density, ref_pos=660.0,
        values=[side_burden, side_burden],
        time_values=[0.0, t_final], g=g_vec[2])

    bc_top = momBC.NeumannBC(
        boundary_name="Top", direction=2,
        density=0.0, ref_pos=0.0,
        values=[over_burden, over_burden],
        time_values=[0.0, t_final], g=g_vec[2])

    bc_cavern = momBC.NeumannBC(
        boundary_name="Cavern", direction=2,
        density=gas_density, ref_pos=430.0,
        values=cavern_values,
        time_values=cavern_times, g=g_vec[2])

    handler = momBC.BcHandler(mom_eq)
    handler.add_boundary_condition(bc_west)
    handler.add_boundary_condition(bc_south)
    handler.add_boundary_condition(bc_bottom)
    handler.add_boundary_condition(bc_east)
    handler.add_boundary_condition(bc_north)
    handler.add_boundary_condition(bc_top)
    handler.add_boundary_condition(bc_cavern)
    return handler


def attach_outputs(mom_eq, folder):
    out = sf.SaveFields(mom_eq)
    out.set_output_folder(folder)
    out.add_output_field("u",        "Displacement (m)")
    out.add_output_field("eps_tot",  "Total strain (-)")
    out.add_output_field("sig",      "Stress (Pa)")
    out.add_output_field("p_elems",  "Mean stress (Pa)")
    out.add_output_field("q_elems",  "Von Mises stress (Pa)")
    return out


def main():
    # Grid
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_irregular")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Momentum equation (theta=0.5 -> Crank-Nicolson)
    mom_eq = sf.LinearMomentumMixed(grid, theta=0.5)

    # Linear solver: same as examples/.../Simulation/Run.py
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("gmres")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Elastic-only material first
    mat, salt_density, E0, nu0 = build_elastic_material(mom_eq)
    mom_eq.set_material(mat)

    # Body force
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)

    # Isothermal temperature
    T_kelvin = 298.0
    T0_field = T_kelvin * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Boundary-load magnitudes
    side_burden = 10.0 * ut.MPa
    over_burden = 10.0 * ut.MPa
    p_max       = 10.0 * ut.MPa
    p_min       = 7.0  * ut.MPa
    gas_density = 0.082

    # Operation BCs: 2 cycles between p_max and p_min over 30 days
    # (each 15-day cycle: 2 d ramp down, 11 d hold low, 2 d ramp up).
    tc_op = sf.TimeController(
        dt=2.0, initial_time=0.0, final_time=30.0, time_unit="day")

    cavern_values = [p_max, p_min, p_min, p_max,
                     p_min, p_min, p_max]
    cavern_times = [0.0,
                    2.0  * ut.day,
                    13.0 * ut.day,
                    15.0 * ut.day,
                    17.0 * ut.day,
                    28.0 * ut.day,
                    30.0 * ut.day]

    bc_op = build_bcs(
        mom_eq, t_final=tc_op.t_final,
        salt_density=salt_density, g_vec=g_vec,
        side_burden=side_burden, over_burden=over_burden,
        cavern_values=cavern_values,
        cavern_times=cavern_times,
        gas_density=gas_density)
    mom_eq.set_boundary_conditions(bc_op)

    # One-shot elastic solve -> populate mom_eq.sig with the initial stress.
    mom_eq.bc.update_dirichlet(0.0)
    mom_eq.bc.update_neumann(0.0)
    mom_eq.bc.update_cavern_bcs(sf.CavernBC.CavernHandler())
    mom_eq.solve_elastic_response()
    eps_tot = mom_eq.compute_total_strain()
    mom_eq.compute_elastic_stress(eps_tot)

    # Build MD, init zeta from elastic stress, then attach
    md = build_munson_dawson(mom_eq, SCENARIO, E0, nu0)
    initialize_md_steady_state(md, mom_eq, T_kelvin=T_kelvin)
    mom_eq.mat.add_to_non_elastic(md)

    out_folder = os.path.join("output", "munson_dawson_example")
    if MPI.COMM_WORLD.rank == 0:
        print(f"[operation] output -> {out_folder}")
    out_op = attach_outputs(mom_eq, out_folder)

    sim_op = sf.Simulator_M(mom_eq, tc_op, [out_op], compute_elastic_response=False)
    sim_op.run()


if __name__ == "__main__":
    main()
