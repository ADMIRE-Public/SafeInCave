import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "safeincave"))
from InputFileAssistant import BuildInputFile
from Simulator import Simulator

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
ifa.set_krylov_solver(method="cg", preconditioner="ilu", rel_tol=1e-12)
# ifa.set_direct_solver(method="petsc")

# Create simulation_settings section
ifa.set_equilibrium_stage(active=False, dt=0.05*hour, tol=1e-9)
ifa.set_operation_stage(active=True, dt=0.05*hour, n_skip=1)


# Create body_forces section
ifa.section_body_forces(density=0.0, direction=2)

# Create time_settings section
time_list = [0*hour,  1*hour]
ifa.section_time(time_list, theta=0.0)

# # Create boundary_conditions section
# ifa.section_boundary_conditions()

# Add Dirichlet boundary conditions
ifa.add_dirichlet(name="WEST", values=list(np.zeros(len(time_list))), component=0)
ifa.add_dirichlet(name="SOUTH", values=list(np.zeros(len(time_list))), component=1)
ifa.add_dirichlet(name="BOTTOM", values=list(np.zeros(len(time_list))), component=2)

# Add Neumann boundary condition
ifa.add_neumann(name="EAST", values=[5*MPa, 5*MPa])
ifa.add_neumann(name="NORTH", values=[5*MPa, 5*MPa])
ifa.add_neumann(name="TOP", values=[8*MPa, 8*MPa])


# Define constitutive model
# Add elastic elements
ifa.add_elastic_element(name="Spring_0", E=[8*GPa, 10*GPa], nu=[0.2, 0.3], active=True, equilibrium=True)

# Add viscoelastic elements
ifa.add_viscoelastic_element(name="KelvinVoigt_0", E=[8*GPa, 5*GPa], nu=[0.35, 0.28], eta=[105e11, 38e11], active=True, equilibrium=False)


# Save input_file.json
ifa.save_input_file("input_file.json")


# Build simulator
sim = Simulator(ifa.input_file)

# Run simulation
sim.run()







