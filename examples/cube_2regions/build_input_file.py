import os
import json
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "safeincave"))
from Grid import GridHandlerGMSH
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
ifa.section_input_grid(path_to_grid, "geom")

# Create output section
ifa.section_output(os.path.join("output", "case_0"))

# Create solver settings section
solver_settings = {
    "type": "KrylovSolver",
    "method": "cg",
    "preconditioner": "petsc_amg",
    "relative_tolerance": 1e-12,
}
ifa.section_solver(solver_settings)

# Create simulation_settings section
ifa.section_simulation(
	simulation_settings = {
		"equilibrium": {
		"active": False,
		"dt_max": 0.5*hour,
		"time_tol": 1e-4
	},
		"operation": {
		"active": True,
		"dt_max": 0.005*hour,
		"n_skip": 1
	}
}
)

# Create body_forces section
salt_density = 2000
ifa.section_body_forces(value=salt_density, direction=2)

# Create time_settings section
time_list = [0*hour,  1*hour]
ifa.section_time(time_list, theta=0.0)

# # Create boundary_conditions section
ifa.section_boundary_conditions()

# Add Dirichlet boundary conditions
ifa.add_boundary_condition(
	boundary_name = "WEST",
	bc_data = {
		"type": "dirichlet",
		"component": 0,
		"values": list(np.zeros(len(time_list)))
	}
)
ifa.add_boundary_condition(
	boundary_name = "SOUTH",
	bc_data = {
		"type": "dirichlet",
		"component": 1,
		"values": list(np.zeros(len(time_list)))
	}
)
ifa.add_boundary_condition(
	boundary_name = "BOTTOM",
	bc_data = {
		"type": "dirichlet",
		"component": 2,
		"values": list(np.zeros(len(time_list)))
	}
)

# Add Neumann boundary condition
ifa.add_boundary_condition(
	boundary_name = "EAST",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": 0*salt_density,
		"reference_position": 1.0,
		"values": [5*MPa, 5*MPa]
	}
)
ifa.add_boundary_condition(
	boundary_name = "NORTH",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": 0*salt_density,
		"reference_position": 1.0,
		"values": [5*MPa, 5*MPa]
	}
)
ifa.add_boundary_condition(
	boundary_name = "TOP",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": 0.0,
		"reference_position": 1.0,
		"values": [8*MPa, 8*MPa]
	}
)

index_A = []
index_B = []

# Sweep over the grid regions and elements
for cell in do.cells(ifa.grid.mesh):
	region_marker = ifa.grid.subdomains[cell]
	if region_marker == ifa.grid.get_subdomain_tags("OMEGA_A"):
		index_A.append(cell.index())
	elif region_marker == ifa.grid.get_subdomain_tags("OMEGA_B"):
		index_B.append(cell.index())
	else:
		raise Exception("Subdomain tag is invalid. Check your mesh file.")

# Assign material properties
ifa.section_constitutive_model()

# Add elastic properties
E = np.zeros(ifa.n_elems)
E[index_A] = 8*GPa
E[index_B] = 10*GPa

nu = np.zeros(ifa.n_elems)
nu[index_A] = 0.2
nu[index_B] = 0.3

ifa.add_elastic_element(	
	element_name = "Spring_0", 
	element_parameters = {
		"type": "Spring",
		"active": True,
		"parameters": {
			"E": list(E),
			"nu": list(nu)
		}
	}
)

# Add viscoelastic properties
E[index_A] = 8*GPa
E[index_B] = 5*GPa

nu[index_A] = 0.35
nu[index_B] = 0.28

eta = np.zeros(ifa.n_elems)
eta[index_A] = 105e11
eta[index_B] = 38e11

# Add viscoelastic properties
ifa.add_viscoelastic_element( 	
	element_name = "KelvinVoigt_0", 
	element_parameters = {
		"type": "KelvinVoigt",
		"active": True,
		"parameters": {
			"E": 	list(E),
			"nu": 	list(nu),
			"eta": 	list(eta)
		}
	}
)

# Save input_file.json
ifa.save_input_file("input_file.json")


