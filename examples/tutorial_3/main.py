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

def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	salt_thickness = float(data[11][len("salt_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall

def create_input_file():
	# Initialize input file object
	ifa = BuildInputFile()

	# Create input_grid section
	path_to_grid = os.path.join("..", "..", "grids", "cavern_overburden")
	ifa.set_input_grid(path_to_grid, "geom")

	# Create output section
	ifa.set_output_folder(os.path.join("output", "case_0"))

	# Create solver settings section
	# ifa.set_krylov_solver(method="cg", preconditioner="petsc_amg", rel_tol=1e-12)
	ifa.set_krylov_solver(method="cg", preconditioner="ilu", rel_tol=1e-12)
	# ifa.set_direct_solver(method="petsc")

	# Create simulation_settings section
	ifa.set_equilibrium_stage(active=True, dt=0.25*day, tol=1e-4)
	ifa.set_operation_stage(active=True, dt=2*hour, n_skip=5)

	# Define densities
	salt_density = 2200
	ovb_density = 2800
	gas_density = 10


	# Create body_forces section
	ifa.section_body_forces(density=[salt_density, ovb_density], direction=2)

	# Create time_settings section
	time_list = [0*day,  2*day,  6*day, 8*day, 10*day]
	ifa.section_time(time_list, theta=0.5)

	# Create boundary_conditions section

	# Add Dirichlet boundary conditions
	ifa.add_dirichlet(name="West_salt", values=list(np.zeros(len(time_list))), component=0)
	ifa.add_dirichlet(name="West_ovb", values=list(np.zeros(len(time_list))), component=0)
	ifa.add_dirichlet(name="South_salt", values=list(np.zeros(len(time_list))), component=1)
	ifa.add_dirichlet(name="South_ovb", values=list(np.zeros(len(time_list))), component=1)
	ifa.add_dirichlet(name="Bottom", values=list(np.zeros(len(time_list))), component=2)
	ifa.add_dirichlet(name="East_salt", values=list(np.zeros(len(time_list))), component=0)
	ifa.add_dirichlet(name="East_ovb", values=list(np.zeros(len(time_list))), component=0)
	ifa.add_dirichlet(name="North_salt", values=list(np.zeros(len(time_list))), component=1)
	ifa.add_dirichlet(name="North_ovb", values=list(np.zeros(len(time_list))), component=1)

	# Extract geometry dimensions
	Lx = ifa.grid.Lx
	Ly = ifa.grid.Ly
	Lz = ifa.grid.Lz
	z_surface = 0.0

	g = 9.81
	ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(path_to_grid)
	cavern_roof = ovb_thickness + hanging_wall
	p_roof = 0 + salt_density*g*hanging_wall + ovb_density*g*ovb_thickness
	print(p_roof/MPa, 0.2*p_roof/MPa, 0.8*p_roof/MPa)

	# Pressure at the top of the salt layer (bottom of overburden)
	p_top = ovb_density*g*ovb_thickness


	# Add Neumann boundary condition
	ifa.add_neumann(name="Cavern", values=[0.8*p_roof, 0.2*p_roof, 0.2*p_roof, 0.8*p_roof, 0.8*p_roof], direction=2, density=gas_density, reference_position=cavern_roof)
	ifa.add_neumann(name="Top", values=[0*MPa, 0*MPa, 0*MPa, 0*MPa, 0*MPa], direction=2, density=0.0, reference_position=z_surface)

	# Define constitutive model
	# Add elastic elements
	ifa.add_elastic_element(
							name="Spring_0",
							E=[102*GPa, 180*GPa],
							nu=0.3,
							active=True,
							equilibrium=True
	)

	# Add viscoelastic elements
	ifa.add_viscoelastic_element(
							name="KelvinVoigt_0",
							E=10*GPa,
							nu=0.32,
							eta=105e11,
							active=False,
							equilibrium=False
	)

	# Add inelastic elements
	ifa.add_dislocation_creep_element(
							name="disCreep",
							A=[1.9e-20, 0.0],
							n=3.0,
							Q=51600,
							T=298,
							active=True,
							equilibrium=True
	)

	ifa.add_desai_element(
							name="desai",
							mu_1=[5.3665857009859815e-11, 0.0],
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
							active=False,
							equilibrium=False
	)


	# Save input_file.json
	ifa.save_input_file("input_file.json")

	return ifa.input_file


def main():
	# Create input file
	input_file = create_input_file()

	# Build simulator
	sim = Simulator(input_file)

	# Run simulation
	sim.run()

if __name__ == '__main__':
	main()



