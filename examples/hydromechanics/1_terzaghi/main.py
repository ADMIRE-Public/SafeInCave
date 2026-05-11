import safeincave as sf
import safeincave.MassBC as massBC
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import dolfinx as do
from safeincave.Utils import MPa, kPa, read_json
import os
import sys

class Props():
    def __init__(self):
        props = read_json("props.json")
        fluid_name = list(props.keys())[0]
        self.mu = props[fluid_name]["Viscosity"]
        self.c_f = props[fluid_name]["Compressibility"]
        self.rho_f = props[fluid_name]["Density"]

        solid_name = list(props.keys())[1]
        self.c_s = props[solid_name]["Compressibility"]
        self.permeability = props[solid_name]["Permeability"]
        self.phi = props[solid_name]["Porosity"]
        self.nu = props[solid_name]["PoissonsRatio"]
        self.G = props[solid_name]["ShearModulus"]
        self.rho_s = props[solid_name]["Density"]
        self.E = 2*self.G*(1 + self.nu)


class Permeability(sf.PermeabilityBase):
    def __init__(self, n_elems: int, p: do.fem.function.Function, props: Props):
        super().__init__(n_elems)
        self.props = props
        self.p = p
        mesh = p.function_space.mesh
        self.DG0_1 = do.fem.functionspace(mesh, ("DG", 0))

    def compute(self) -> to.Tensor:
        # p_elems = sf.Utils.project(self.p, self.DG0_1)
        # print(self.p.x.array.max())
        # p_elems = to.from_numpy(p_elems.x.array)
        # k = to.where(p_elems < 0.5, 1e-13, 1e-12)
        # k = 1e-12 * to.ones(self.n_elems, dtype=to.float64) + p_elems
        k = (self.props.permeability/self.props.mu) * to.ones(self.n_elems, dtype=to.float64)
        # k = (1 + p_elems**2) * 1e-12
        return k


def main():
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cube")
    grid = sf.GridHandlerGMSH("geom", grid_path)

    # Read poromechanical properties
    props = Props()
    print(props.permeability/props.mu)

    # Define output folder
    output_folder = os.path.join("output", "case_g_stab_1")
    # output_folder = os.path.join("output", "case_g_nonstab")


    # Time settings for equilibrium stage
    t_eq = sf.TimeController(dt=0.5, initial_time=0.0, final_time=5, time_unit="second")


    # Build material properties
    mat = sf.Material(grid.n_elems)

    # Set material density
    rho_m = props.phi*props.rho_f + (1 - props.phi)*props.rho_s
    rho = rho_m*to.ones(grid.n_elems, dtype=to.float64)
    mat.set_solid_density(rho)

    # Set fluid density
    rho_f = props.rho_f*to.ones(grid.n_elems, dtype=to.float64)
    mat.set_fluid_density(rho_f)

    # Set solid grain compressibility
    cs = props.c_s*to.ones(grid.n_elems, dtype=to.float64)
    mat.set_cs(cs)

    # Set fluid compressibility
    cf = props.c_f*to.ones(grid.n_elems, dtype=to.float64)
    mat.set_cf(cf)

    # Set porosity
    phi = props.phi*to.ones(grid.n_elems, dtype=to.float64)
    mat.set_porosity(phi)

    # Set gravitational vector
    mat.set_gravity([0.0, 0.0, -9.81])

    # Set elastic strain model
    E = props.E*to.ones(grid.n_elems, dtype=to.float64)
    nu = props.nu*to.ones(grid.n_elems, dtype=to.float64)
    spring = sf.Spring(E, nu, "spring")
    mat.add_to_elastic(spring)

    # Create poroelastic strain model
    K = E/(3*(1-2*nu))
    G = 3*K*E/(9*K - E)
    biot = 1 - cs*K
    pore_p = sf.Poroelastic(alpha=biot, K=K, name="poroelastic")
    mat.add_to_poroelastic(pore_p)



    ######### Define Mass equation #########
    mass_eq = sf.MassDeformablePorousMedia(grid, is_linear=True, fixed_stress=True)

    # Define solver
    solver_mass = PETSc.KSP().create(grid.mesh.comm)
    solver_mass.setType("gmres")
    solver_mass.getPC().setType("asm")
    solver_mass.setTolerances(rtol=1e-12, max_it=100)
    mass_eq.set_solver(solver_mass)

    # Set permeability
    kappa_model = Permeability(grid.n_elems, mass_eq.P, props)
    mat.set_permeability_model(kappa_model)

    # Set material properties to mass equation
    mass_eq.set_material(mat)

    # Define boundary conditions for mass equation
    bc_handler = massBC.BcHandler(mass_eq)
    mass_eq.set_boundary_conditions(bc_handler)

    # Set initial temperature field
    P0_field = 0.0*to.ones(mass_eq.n_nodes, dtype=to.float64)
    mass_eq.set_initial_P(P0_field)
    mass_eq.update_P_old()


    ######### Define Momentum equation #########
    # mom_eq = sf.LinearMomentumMixed(grid, 0.5)
    mom_eq = sf.LinearMomentum(grid, 0.5)

    # Define solver
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("gmres")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-12, max_it=100)
    mom_eq.set_solver(mom_solver)

    # Set initial temperature field
    T0_field = 293*to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)

    # Set initial pore pressure field
    mom_eq.set_P0(mass_eq.P)
    mom_eq.set_P(mass_eq.P)

    # Set material
    mom_eq.set_material(mat)

    # Boundary conditions
    bc_west = momBC.DirichletBC(boundary_name = "WEST", 
                                component = 0,
                                values = [0.0, 0.0],
                                time_values = [t_eq.t_initial, t_eq.t_final])

    bc_east = momBC.DirichletBC(boundary_name = "EAST", 
                                component = 0,
                                values = [0.0, 0.0],
                                time_values = [t_eq.t_initial, t_eq.t_final])

    bc_south = momBC.DirichletBC(boundary_name = "SOUTH", 
                                component = 1,
                                values = [0.0, 0.0],
                                time_values = [t_eq.t_initial, t_eq.t_final])

    bc_north = momBC.DirichletBC(boundary_name = "NORTH", 
                                component = 1,
                                values = [0.0, 0.0],
                                time_values = [t_eq.t_initial, t_eq.t_final])

    bc_bottom = momBC.DirichletBC(boundary_name = "BOTTOM", 
                                component = 2,
                                values = [0.0, 0.0],
                                time_values = [t_eq.t_initial, t_eq.t_final])

    sigma_z = 10.0*kPa
    bc_top = momBC.NeumannBC(boundary_name = "TOP",
                                direction = 0,
                                density = 0,
                                ref_pos = 0,
                                values = [sigma_z, sigma_z],
                                time_values = [t_eq.t_initial, t_eq.t_final])

    bc_handler = momBC.BcHandler(mom_eq)
    bc_handler.add_boundary_condition(bc_west)
    bc_handler.add_boundary_condition(bc_south)
    bc_handler.add_boundary_condition(bc_bottom)
    bc_handler.add_boundary_condition(bc_east)
    bc_handler.add_boundary_condition(bc_north)
    bc_handler.add_boundary_condition(bc_top)

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_handler)

    # Bind displacement field to mass equation for Biot coupling
    mass_eq.update_u_old(mom_eq.u)

    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_folder_eq = os.path.join(output_folder, "equilibrium")
    output_mom.set_output_folder(output_folder_eq)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")

    output_mass = sf.SaveFields(mass_eq)
    output_mass.set_output_folder(output_folder_eq)
    output_mass.add_output_field("P", "Pressure (Pa)")

    outputs = [output_mom, output_mass]

    # Define simulator
    sim = sf.Simulator_HM(mom_eq, mass_eq, t_eq, outputs, compute_elastic_response=True)
    sim.run()



    # ###############################################################
    # #################### RUN OPERATIONAL STAGE ####################
    # ###############################################################


    # # Time settings for equilibrium stage
    # t_op = sf.TimeController(dt=0.05, initial_time=0.0, final_time=5.0, time_unit="second")

    # # Define boundary conditions for mass equation
    # bc_handler = massBC.BcHandler(mass_eq)

    # # Define boundary conditions for mass equation
    # bc_top = massBC.DirichletBC(boundary_name = "TOP", 
    #                         values = [0.0, 0.0],
    #                         time_values = [t_op.t_initial, t_op.t_final])

    # bc_handler = massBC.BcHandler(mass_eq)
    # bc_handler.add_boundary_condition(bc_top)
    # mass_eq.set_boundary_conditions(bc_handler)

    # # Create output handlers
    # output_mom = sf.SaveFields(mom_eq)
    # output_folder_op = os.path.join(output_folder, "operation")
    # output_mom.set_output_folder(output_folder_op)
    # output_mom.add_output_field("u", "Displacement (m)")
    # output_mom.add_output_field("sig", "Stress (Pa)")
    # output_mom.add_output_field("p_elems", "Mean stress (Pa)")

    # output_mass = sf.SaveFields(mass_eq)
    # output_mass.set_output_folder(output_folder_op)
    # output_mass.add_output_field("P", "Pressure (Pa)")

    # outputs = [output_mom, output_mass]

    # # Define simulator
    # sim = sf.Simulator_HM_twoway(mom_eq, mass_eq, t_op, outputs, False)
    # sim.run()



if __name__ == '__main__':
	main()