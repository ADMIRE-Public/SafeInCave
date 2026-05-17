"""
Minimal working example: cube under triaxial loading with the Modified
Cam-Clay model.

Model stack:
    Spring (linear elastic) + ModifiedCamClayViscoplastic

Geometry: grids/cube (symmetry-quadrant of a unit-ish cube).
Schedule:
    0 - 1 h : isotropic consolidation, sigma_1 = sigma_3 = 10 MPa
    1 - 3 h : load   axial sigma_1 -> 25 MPa, lateral sigma_3 held at 10 MPa
    3 - 5 h : unload axial sigma_1 -> 10 MPa
    5 - 7 h : reload axial sigma_1 -> 25 MPa
    7 - 9 h : unload axial sigma_1 -> 10 MPa
The two load/unload cycles exercise the cyclic-hardening update of p_c.

Material parameters: Brine 100 bar baseline from Saeed's notebook
(examples/cam_clay/Calibration_Clay_Rich.ipynb, ``table_guess``).

Run: python main.py  -> outputs in ./output/case_0/<formulation>/
"""

import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
import dolfinx as do
import os
import torch as to


class LinearMomentumMod(sf.LinearMomentum):
    def __init__(self, grid, theta):
        super().__init__(grid, theta)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.pc = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[0].eps_ne_k)
        self.Fvp.x.array[:] = self.mat.elems_ne[0].Fvp
        self.pc.x.array[:] = self.mat.elems_ne[0].pc


class LinearMomentumMixedMod(sf.LinearMomentumMixed):
    def __init__(self, grid, theta, stab_scaling=1.0):
        super().__init__(grid, theta, stab_scaling)
        self.Fvp = do.fem.Function(self.DG0_1)
        self.pc = do.fem.Function(self.DG0_1)
        self.eps_vp = do.fem.Function(self.DG0_3x3)

    def run_after_solve(self):
        self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[0].eps_ne_k)
        self.Fvp.x.array[:] = self.mat.elems_ne[0].Fvp
        self.pc.x.array[:] = self.mat.elems_ne[0].pc


def run(formulation):
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cube")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Output folder
    output_folder = os.path.join("output", "case_0", f"{formulation}")

    # Time settings (1h consolidation + 2 load/unload cycles of 4h each)
    unit = "hour"
    dt = 0.1
    t_final = 9.0
    t_control = sf.TimeController(dt=dt, initial_time=0.0,
                                  final_time=t_final, time_unit=unit)

    # Momentum equation
    theta = 0.5
    if formulation == "P1":
        mom_eq = LinearMomentumMod(grid, theta=theta)
    elif formulation == "P1P1":
        mom_eq = LinearMomentumMixedMod(grid, theta=theta, stab_scaling=0.0)
    elif formulation == "P1P1_Stab":
        mom_eq = LinearMomentumMixedMod(grid, theta=theta, stab_scaling=1.0)

    # Solver
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("gmres")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-10, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Material
    mat = sf.Material(mom_eq.n_elems)

    rho = 2000.0 * to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)

    # Elastic skeleton (clay-rich -> softer than salt)
    E = 2.0 * ut.GPa * to.ones(mom_eq.n_elems)
    nu = 0.30 * to.ones(mom_eq.n_elems)
    spring = sf.Spring(E, nu, "spring")

    # Modified Cam-Clay (Brine 100 bar — Saeed's table_guess)
    phi0 = 0.32
    e0_val = phi0 / (1.0 - phi0)
    mcc = sf.ModifiedCamClayViscoplastic(
        M     = 0.702   * to.ones(mom_eq.n_elems),
        lam   = 0.00285 * to.ones(mom_eq.n_elems),
        kap   = 0.00055 * to.ones(mom_eq.n_elems),
        theta = 0.0045  * to.ones(mom_eq.n_elems),
        pc0   = 13.0 * ut.MPa * to.ones(mom_eq.n_elems),
        e0    = e0_val  * to.ones(mom_eq.n_elems),
        eta_v = 1e11    * to.ones(mom_eq.n_elems),
        n_rate= 2.0     * to.ones(mom_eq.n_elems),
        name  = "cam_clay",
    )

    mat.add_to_elastic(spring)
    mat.add_to_non_elastic(mcc)
    mom_eq.set_material(mat)

    # No body force in a triaxial test
    g_vec = [0.0, 0.0, 0.0]
    mom_eq.build_body_force(g_vec)

    # Isothermal
    T0_field = 293.0 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Boundary conditions — symmetry quadrant of a cube triaxial.
    sigma_3       = 10.0 * ut.MPa   # confining (held constant)
    sigma_1_low   = 10.0 * ut.MPa
    sigma_1_peak  = 25.0 * ut.MPa

    # Schedule waypoints (1h consolidation + 2 load/unload cycles)
    times_h    = [0.0, 1.0, 3.0, 5.0, 7.0, 9.0]
    sigma_1_vs = [sigma_1_low, sigma_1_low,
                  sigma_1_peak, sigma_1_low,
                  sigma_1_peak, sigma_1_low]
    times_s = [t * ut.hour for t in times_h]

    bc_west = momBC.DirichletBC("WEST", 0, [0.0, 0.0],
                                [0.0, t_control.t_final])
    bc_south = momBC.DirichletBC("SOUTH", 1, [0.0, 0.0],
                                 [0.0, t_control.t_final])
    bc_bottom = momBC.DirichletBC("BOTTOM", 2, [0.0, 0.0],
                                  [0.0, t_control.t_final])
    bc_east = momBC.NeumannBC("EAST", direction=2, density=0.0, ref_pos=0.0,
                              values=[sigma_3, sigma_3],
                              time_values=[0.0, t_control.t_final],
                              g=g_vec[2])
    bc_north = momBC.NeumannBC("NORTH", direction=2, density=0.0, ref_pos=0.0,
                               values=[sigma_3, sigma_3],
                               time_values=[0.0, t_control.t_final],
                               g=g_vec[2])
    bc_top = momBC.NeumannBC("TOP", direction=2, density=0.0, ref_pos=0.0,
                             values=sigma_1_vs,
                             time_values=times_s,
                             g=g_vec[2])

    bc_handler = momBC.BcHandler(mom_eq)
    bc_handler.add_boundary_condition(bc_west)
    bc_handler.add_boundary_condition(bc_south)
    bc_handler.add_boundary_condition(bc_bottom)
    bc_handler.add_boundary_condition(bc_east)
    bc_handler.add_boundary_condition(bc_north)
    bc_handler.add_boundary_condition(bc_top)
    mom_eq.set_boundary_conditions(bc_handler)

    # Outputs
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u",       "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("sig",     "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    output_mom.add_output_field("eps_vp",  "Viscoplastic strain (-)")
    output_mom.add_output_field("Fvp",     "Yield function F (Pa^2)")
    output_mom.add_output_field("pc",      "Preconsolidation pressure (Pa)")
    outputs = [output_mom]

    sim = sf.Simulator_M(mom_eq, t_control, outputs,
                         compute_elastic_response=True)
    sim.run()


def main():
    run("P1P1")
    # run("P1")
    # run("P1P1_Stab")


if __name__ == "__main__":
    main()
