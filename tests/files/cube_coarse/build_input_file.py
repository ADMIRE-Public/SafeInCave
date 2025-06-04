import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "..", "safeincave"))
from Grid import GridHandlerGMSH
from InputFileAssistant import BuildInputFile

# Useful units
hour = 60*60
day = 24*hour
MPa = 1e6

# Initialize input file object
bif = BuildInputFile()

# Create input_grid section
# path_to_grid = os.path.join("files", "cube_coarse")
bif.set_input_grid(".", "geom")

# Create output section
bif.set_output_folder(os.path.join("output", "case_0"))

# Create solver settings section
bif.set_solver(solver_type="cg", solver_PC="gamg", rtol=1e-12, maxite=100)

# Create simulation_settings section
bif.set_equilibrium_stage(active=False, dt=0.5*hour, tol=1e-4, ite_max=50)
bif.set_operation_stage(active=True, dt=0.1*hour, n_skip=1)

# Create body_forces section
salt_density = 2000
bif.section_body_forces(density=salt_density, direction=2)

# Create time_settings section
time_list = [0*hour,  2*hour,  14*hour, 16*hour, 24*hour]
bif.section_time(time_list, theta=0.0)

# # Create boundary_conditions section
bif.section_boundary_conditions()

# Add Dirichlet boundary conditions
bif.add_dirichlet(name="WEST", values=list(np.zeros(len(time_list))), component=0)
bif.add_dirichlet(name="SOUTH", values=list(np.zeros(len(time_list))), component=1)
bif.add_dirichlet(name="BOTTOM", values=list(np.zeros(len(time_list))), component=2)

# Add Neumann boundary condition
bif.add_neumann(name="EAST", values=[4*MPa, 4*MPa, 4*MPa, 4*MPa, 4*MPa], direction=2, density=0.0, reference_position=1.0)
bif.add_neumann(name="NORTH", values=[4*MPa, 4*MPa, 4*MPa, 4*MPa, 4*MPa], direction=2, density=0.0, reference_position=1.0)
bif.add_neumann(name="TOP", values=[4.1*MPa, 12*MPa, 12*MPa, 6*MPa, 6*MPa], direction=2, density=0.0, reference_position=1.0)


# Build constitutive model
bif.add_elastic_element(
						name="Spring_0",
						E=102e9,
						nu=0.3,
						active=True,
						equilibrium=True
)

bif.add_viscoelastic_element(
						name="KelvinVoigt_0",
						E=10e9,
						nu=0.32,
						eta=105e11,
						active=True,
						equilibrium=False
)

bif.add_desai_element(
						name="desai",
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
						active=True,
						equilibrium=False
)

bif.add_dislocation_creep_element(
						name="disCreep",
						A=1.9e-20,
						n=3.0,
						Q=51600,
						T=298,
						active=True,
						equilibrium=False
)


# Save input_file.json
bif.save_input_file("input_file.json")

