import os
import sys
import numpy as np
sys.path.append(os.path.join("..", "..", "safeincave"))
from Grid import GridHandlerGMSH
from InputFileAssistant import BuildInputFile

# Useful units
hour = 60*60
day = 24*hour
MPa = 1e6

# Initialize input file object
bif = BuildInputFile()

# Create input_grid section
path_to_grid = os.path.join("..", "..", "grids", "cavern_irregular")
bif.section_input_grid(path_to_grid, "geom")

print(bif.grid.get_boundary_names())
print(bif.grid.get_subdomain_names())

# Extract geometry dimensions
Lx = bif.grid.Lx
Ly = bif.grid.Ly
Lz = bif.grid.Lz
cavern_roof = 430

# Create output section
bif.section_output(os.path.join("output", "case_0"))

# Create solver settings section
solver_settings = {
    "type": "KrylovSolver",
    "method": "cg",
    "preconditioner": "petsc_amg",
    "relative_tolerance": 1e-12,
}
bif.section_solver(solver_settings)

# Create simulation_settings section
bif.section_simulation(
	simulation_settings = {
		"equilibrium": {
			"active": True,
			"dt_max": 0.5*hour,
			"time_tol": 1e-4
		},
		"operation": {
			"active": True,
			"dt_max": 0.05*hour,
			"n_skip": 4
		}
	}
)

# Create body_forces section
salt_density = 2000
bif.section_body_forces(value=salt_density, direction=2)

# Create time_settings section
time_list = [0*hour,  2*hour,  14*hour, 16*hour, 24*hour]
bif.section_time(time_list, theta=0.0)

# Create boundary_conditions section
bif.section_boundary_conditions()

# Add Dirichlet boundary conditions
bif.add_boundary_condition(
	boundary_name = "West",
	bc_data = {
		"type": "dirichlet",
		"component": 0,
		"values": list(np.zeros(len(time_list)))
	}
)
bif.add_boundary_condition(
	boundary_name = "South",
	bc_data = {
		"type": "dirichlet",
		"component": 1,
		"values": list(np.zeros(len(time_list)))
	}
)
bif.add_boundary_condition(
	boundary_name = "Bottom",
	bc_data = {
		"type": "dirichlet",
		"component": 2,
		"values": list(np.zeros(len(time_list)))
	}
)

# Add Neumann boundary condition
bif.add_boundary_condition(
	boundary_name = "East",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": salt_density,
		"reference_position": Lz,
		"values": [10*MPa, 10*MPa, 10*MPa, 10*MPa, 10*MPa]
	}
)

bif.add_boundary_condition(
	boundary_name = "North",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": salt_density,
		"reference_position": Lz,
		"values": [10*MPa, 10*MPa, 10*MPa, 10*MPa, 10*MPa]
	}
)

bif.add_boundary_condition(
	boundary_name = "Top",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": salt_density,
		"reference_position": Lz,
		"values": [10*MPa, 10*MPa, 10*MPa, 10*MPa, 10*MPa]
	}
)

h2_density = 10
bif.add_boundary_condition(
	boundary_name = "Cavern",
	bc_data = {
		"type": "neumann",
		"direction": 2,
		"density": h2_density,
		"reference_position": cavern_roof,
		"values": [10*MPa, 7*MPa, 7*MPa, 10*MPa, 10*MPa]
	}
)

# Assign material properties
bif.section_constitutive_model()

# Add elastic properties
bif.add_elastic_element(	
	element_name = "Spring0", 
	element_parameters = {
		"type": "Spring",
		"active": True,
		"parameters": {
			"E":  list(102e9*np.ones(bif.n_elems)),
			"nu": list(0.3*np.ones(bif.n_elems))
		}
	}
)

# Add viscoelastic properties
bif.add_viscoelastic_element( 	
	element_name = "KelvinVoigt1", 
	element_parameters = {
		"type": "KelvinVoigt",
		"active": True,
		"parameters": {
			"E":   list(10e9*np.ones(bif.n_elems)),
			"nu":  list(0.32*np.ones(bif.n_elems)),
			"eta": list(105e11*np.ones(bif.n_elems))
		}
	}
)

# Add viscoplastic parameters
bif.add_inelastic_element(	
	element_name = "ViscPlastDesai", 
	element_parameters = {
		"type": "ViscoplasticDesai",
		"active": True,
		"parameters": {
			"mu_1": 	list(5.3665857009859815e-11*np.ones(bif.n_elems)),
			"N_1": 		list(3.1*np.ones(bif.n_elems)),
			"n": 		list(3.0*np.ones(bif.n_elems)),
			"a_1":		list(1.965018496922832e-05*np.ones(bif.n_elems)),
			"eta": 		list(0.8275682807874163*np.ones(bif.n_elems)),
			"beta_1": 	list(0.0048*np.ones(bif.n_elems)),
			"beta": 	list(0.995*np.ones(bif.n_elems)),
			"m": 		list(-0.5*np.ones(bif.n_elems)),
			"gamma": 	list(0.095*np.ones(bif.n_elems)),
			"alpha_0": 	list(0.0022*np.ones(bif.n_elems)),
			"k_v": 		list(0.0*np.ones(bif.n_elems)),
			"sigma_t": 	list(5.0*np.ones(bif.n_elems))
		}
	}
)

# Add dislocation creep parameters
bif.add_inelastic_element(	
	element_name = "DisCreep", 
	element_parameters = {
		"type": "DislocationCreep",
		"active": True,
		"parameters": {
			"A": list(1.9e-20*np.ones(bif.n_elems)),
			"n": list(3.0*np.ones(bif.n_elems)),
			"T": list(298*np.ones(bif.n_elems)),
			"Q": list(51600*np.ones(bif.n_elems)),
			"R": list(8.32*np.ones(bif.n_elems))
		}
	}
)

# Save input_file.json
bif.save_input_file("input_file.json")

