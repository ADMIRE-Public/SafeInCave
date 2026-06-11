# Example 4: Single Cavern With Operation Stage

## Introduction

This example extends the single-cavern mechanical workflow by separating the calculation into an equilibrium stage and a later operation stage. It uses salt creep and viscoelasticity to compute stress and deformation around the cavern.

## Problem description

This section walks through `examples/mechanics/4_cavern/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, PETSc linear solvers, field sampling utilities. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
from safeincave.Utils import day, GPa, create_field_elems
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import os
import sys
```

### Block 2: Function `get_geometry_parameters`

This helper reads geometric constants directly from the mesh-generation file. The example uses those dimensions to derive physically meaningful pressure references, instead of hard-coding cavern roof or overburden values in several places.

```Python

def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, hanging_wall
```

### Block 3: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python


def run(formulation):
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_overburden_coarse")
	# grid_path = os.path.join("..", "..", "grids", "cavern_overburden")
	grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 4: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder = os.path.join("output", "case_0", f"{formulation}")
```

### Block 5: Define momentum equation

The mixed momentum formulation is selected when pressure-related stability or volumetric behavior matters. The `stab_scaling` parameter controls the pressure-stabilization contribution, while `theta` sets the time-integration convention used by inelastic updates.

```Python
	# Define momentum equation
	theta = 0.0
	if formulation == "P1":
		mom_eq = sf.LinearMomentum(grid, theta=theta)
	elif formulation == "P1P1":
		mom_eq = sf.LinearMomentumMixed(grid, theta=theta, stab_scaling=0.0)
	elif formulation == "P1P1_Stab":
		mom_eq = sf.LinearMomentumMixed(grid, theta=theta, stab_scaling=1.0)
```

### Block 6: Define solver

The mechanical linear system is handed to PETSc. GMRES is used for the momentum formulations because mixed forms, stabilization, and history-dependent tangents can make the assembled system less friendly to simple symmetric solvers.

```Python
	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("gmres")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)
```

### Block 7: Define material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python
	# Define material properties
	mat = sf.Material(mom_eq.n_elems)
```

### Block 8: Extract region indices

The region-index dictionaries come from `GridHandlerGMSH`. They let the script assign arrays element by element, so salt, overburden, or separate cube regions can share one mesh while keeping distinct properties.

```Python
	# Extract region indices
	ind_salt = grid.region_indices["Salt"]
	ind_ovb = grid.region_indices["Overburden"]
```

### Block 9: Set material density

Density is assigned as an element-wise tensor rather than a scalar so each geological region can carry its own body-force contribution. The same vector is also available to material routines that need density-dependent quantities.

```Python
	# Set material density
	salt_density = 2200
	ovb_density = 2800
	gas_density = 0.082
	rho = to.zeros(mom_eq.n_elems, dtype=to.float64)
	rho[ind_salt] = salt_density
	rho[ind_ovb] = ovb_density
	mat.set_density(rho)
```

### Block 10: Constitutive model

The spring is the elastic part of the constitutive model. Here its Young modulus and Poisson ratio are assembled as element-wise arrays, which is why different mesh regions can receive different elastic stiffness.

```Python
	# Constitutive model
	E0 = to.zeros(mom_eq.n_elems)
	E0[ind_salt] = 102*GPa
	E0[ind_ovb] = 180*GPa
	nu0 = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E0, nu0, "spring")
```

### Block 11: Create Kelvin-Voigt viscoelastic element

The Kelvin-Voigt `Viscoelastic` element adds delayed, rate-dependent strain. It is included when the example needs short-term transient deformation beyond purely elastic response.

```Python
	# Create Kelvin-Voigt viscoelastic element
	eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)
	kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")
```

### Block 12: Create creep

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
	# Create creep
	A = to.zeros(mom_eq.n_elems)
	A[ind_salt] = 1.9e-20
	A[ind_ovb] = 0.0
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_0 = sf.DislocationCreep(A, Q, n, "creep")
```

### Block 13: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
	# Create constitutive model
	mat.add_to_elastic(spring_0)
	# mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_0)
```

### Block 14: Set constitutive model

Attaching the material to the momentum equation initializes stiffness-related fields and makes the constitutive mechanisms available to the mechanical solver.

```Python
	# Set constitutive model
	mom_eq.set_material(mat)
```

### Block 15: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)
```

### Block 16: Set initial temperature field

The temperature function is evaluated at element centroids. This is the right representation for mechanical constitutive laws because creep and thermoelastic strain are stored per element.

```Python
	# Set initial temperature field
	def T_field_fun(x,y,z):
		km = 1000
		dTdZ = 27/km
		T_surface = 20 + 273
		return T_surface - dTdZ*z
	T0_field = create_field_elems(grid, T_field_fun)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)
```

### Block 17: Time settings for equilibrium stage

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
	# Time settings for equilibrium stage
	tc_eq = sf.TimeControllerParabolic(n_time_steps=20,
										initial_time=0.0,
										final_time=5,
										time_unit="day")
```

### Block 18: Boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Boundary conditions
	bc_equilibrium = momBC.BcHandler(mom_eq)
```

### Block 19: Apply Dirichlet boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	# Apply Dirichlet boundary conditions
	boundaries = [("West_salt", 0),
					("West_ovb", 0),
					("East_salt", 0),
					("East_ovb", 0),
					("South_salt", 1),
					("South_ovb", 1),
					("North_salt", 1),
					("North_ovb", 1),
					("Bottom", 2)]
	for b_name, component in boundaries:
		bc = momBC.DirichletBC(boundary_name=b_name, component=component, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
		bc_equilibrium.add_boundary_condition(bc)
```

### Block 20: Extract geometry dimensions

The mesh is loaded through `GridHandlerGMSH`, which does more than read coordinates and cells. It also exposes named facet tags for boundary conditions, region tags for material assignment, mesh dimensions, cell volumes, and smoothing operators used when element stresses are reported at nodes or smoothed back to elements.

```Python
	# Extract geometry dimensions
	ovb_thickness, hanging_wall = get_geometry_parameters(grid_path)
```

### Block 21: Calculate lithostatic pressure at the cavern's roof

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Calculate lithostatic pressure at the cavern's roof
	cavern_roof = ovb_thickness + hanging_wall
	p_roof = 0 + salt_density*abs(g)*hanging_wall + ovb_density*abs(g)*ovb_thickness
	print(0.2*p_roof/1e6, 0.8*p_roof/1e6, cavern_roof)
```

### Block 22: Impose loading condition on the cavern walls

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	# Impose loading condition on the cavern walls
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = [0.8*p_roof, 0.8*p_roof],
						time_values = [0*day,  tc_eq.t_final],
						g = g_vec[2])
	bc_equilibrium.add_boundary_condition(bc_cavern)
```

### Block 23: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python

	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_equilibrium)
```

### Block 24: Equilibrium output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Equilibrium output folder
	ouput_folder_equilibrium = os.path.join(output_folder, "equilibrium")
```

### Block 25: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(ouput_folder_equilibrium)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]
```

### Block 26: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
	# Define simulator
	sim = sf.Simulator_M(mom_eq, tc_eq, outputs, compute_elastic_response=True)
	sim.run()
```

### Block 27: Time settings for operation stage

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python






	# Time settings for operation stage
	tc_op = sf.TimeController(dt=2, initial_time=0.0, final_time=240, time_unit="hour")
```

### Block 28: Boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Boundary conditions
	bc_operation = momBC.BcHandler(mom_eq)
```

### Block 29: Dirichlet boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	# Dirichlet boundary conditions
	boundaries = [("West_salt", 0), ("West_ovb", 0), ("East_salt", 0), ("East_ovb", 0),
				  ("South_salt", 1), ("South_ovb", 1), ("North_salt", 1), ("North_ovb", 1),
				  ("Bottom", 2)]
	for b_name, component in boundaries:
		bc = momBC.DirichletBC(boundary_name=b_name, component=component, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
		bc_operation.add_boundary_condition(bc)
```

### Block 30: Impose loading condition on the cavern walls

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	# Impose loading condition on the cavern walls
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = [0.8*p_roof, 0.2*p_roof, 0.2*p_roof, 0.8*p_roof, 0.8*p_roof],
						time_values = [0*day,  2*day,  6*day, 8*day, 10*day],
						g = g_vec[2])
	bc_operation.add_boundary_condition(bc_cavern)
```

### Block 31: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_operation)
```

### Block 32: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder_operation = os.path.join(output_folder, "operation")
```

### Block 33: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder_operation)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]
```

### Block 34: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
	# Define simulator
	sim = sf.Simulator_M(mom_eq, tc_op, outputs, compute_elastic_response=False)
	sim.run()
```

### Block 35: run("P1P1")

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
	run("P1")
	# run("P1P1")
	# run("P1P1_Stab")
```

### Block 36: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python
if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: loads the cavern mesh.
- `safeincave.LinearMomentum` and `safeincave.LinearMomentumMixed`: define alternative mechanical formulations.
- `safeincave.Material`: stores material density and constitutive elements.
- `safeincave.Spring`, `safeincave.Viscoelastic`, and `safeincave.DislocationCreep`: define the salt mechanical model.
- `safeincave.TimeControllerParabolic` and `safeincave.TimeController`: control equilibrium and operation stages.
- `safeincave.MomentumBC.DirichletBC` and `safeincave.MomentumBC.NeumannBC`: apply displacement constraints, overburden loading, and cavern pressure.
- `safeincave.Simulator_M`: runs each mechanical stage.

## Running the Example

From the repository root, run:

```bash
cd examples/mechanics/4_cavern
python main.py
```

The folder also contains helper scripts for alternative runs and post-processing:

```bash
python plot_results.py
python plot_cavern_data.py
python plot_volume_loss.py
```

## Conclusion

This example demonstrates a staged cavern mechanics analysis with equilibrium initialization followed by operation. It can be adapted for pressure histories, convergence comparisons, and long-term creep studies.
