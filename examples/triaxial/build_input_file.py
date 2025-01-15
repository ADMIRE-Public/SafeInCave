import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "safeincave"))
from InputFileAssistant import BuildInputFile
import dolfin as do

# Useful units
hour = 60*60
day = 24*hour
MPa = 1e6
GPa = 1e9

# Initialize input file object
ifa = BuildInputFile()

# Create input_grid section
path_to_grid = os.path.join("..", "..", "grids", "cube")
ifa.set_input_grid(path_to_grid, "geom")

# Create output section
ifa.set_output_folder(os.path.join("output", "case_0"))

# Create solver settings section
ifa.set_krylov_solver(method="cg", preconditioner="petsc_amg", rel_tol=1e-12)
# ifa.set_direct_solver(method="petsc")

# Create simulation_settings section
ifa.set_equilibrium_stage(active=True, dt=0.5*hour, tol=1e-4)
ifa.set_operation_stage(active=True, dt=0.1*hour, n_skip=2)


# Create body_forces section
ifa.section_body_forces(density=0.0, direction=2)

# Create time_settings section
time_list = [0*hour,  2*hour,  14*hour, 16*hour, 24*hour]
ifa.section_time(time_list, theta=0.5)

# # Create boundary_conditions section
# ifa.section_boundary_conditions()

# Add Dirichlet boundary conditions
ifa.add_dirichlet(name="WEST", values=list(np.zeros(len(time_list))), component=0)
ifa.add_dirichlet(name="SOUTH", values=list(np.zeros(len(time_list))), component=1)
ifa.add_dirichlet(name="BOTTOM", values=list(np.zeros(len(time_list))), component=2)

# Add Neumann boundary condition
ifa.add_neumann(name="EAST", values=[4*MPa, 4*MPa, 4*MPa, 4*MPa, 4*MPa])
ifa.add_neumann(name="NORTH", values=[4*MPa, 4*MPa, 4*MPa, 4*MPa, 4*MPa])
ifa.add_neumann(name="TOP", values=[4.1*MPa, 12*MPa, 12*MPa, 6*MPa, 6*MPa])


# Define constitutive model

# Add elastic elements
ifa.add_elastic_element(name="Spring_0", E=102*GPa, nu=0.3, active=True)

# Add viscoelastic elements
ifa.add_viscoelastic_element(name="KelvinVoigt_0", E=10*GPa, nu=0.32, eta=105e11, active=True)

# Add inelastic elements
ifa.add_dislocation_creep_element(name="disCreep", A=1.9e-20, n=3.0, Q=51600, T=298, active=True)

ifa.add_desai_element(	name="desai",
						mu_1=5.3665857009859815e-11,
						N_1=3.1,
						n=3.0,
						a_1=1.965018496922832e-05,
						eta=0.8275682807874163,
						beta_1=0.0048,
						beta=0.995,
						m=-0.5,
						gamma=0.095,
						alpha_0=0.0022,
						sigma_t=5.0,
						active=True )


# Save input_file.json
ifa.save_input_file("input_file.json")


