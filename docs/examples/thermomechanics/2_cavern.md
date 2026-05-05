# Example 2: Thermomechanical Single Cavern

## Introduction

This example models a salt cavern with overburden in two stages. First, it computes a mechanical equilibrium state under geostatic loading and cavern pressure. Then it couples heat diffusion and mechanics during a 240-day operation stage with changing cavern pressure.

## Problem description

This section walks through `examples/thermomechanics/2_cavern/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, thermal boundary conditions, PETSc linear solvers, field sampling utilities. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
import safeincave.Utils as ut
from safeincave.Utils import GPa, MPa, day, hour, create_field_elems, create_field_nodes
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import os
import sys
import torch as to
```

### Block 2: Function `get_geometry_parameters`

This helper reads geometric constants directly from the mesh-generation file. The example uses those dimensions to derive physically meaningful pressure references, instead of hard-coding cavern roof or overburden values in several places.

```Python

def get_geometry_parameters(path_to_grid):
	f = open(os.path.join(path_to_grid, "geom.geo"), "r")
	data = f.readlines()
	ovb_thickness = float(data[10][len("ovb_thickness = "):-2])
	salt_thickness = float(data[11][len("salt_thickness = "):-2])
	hanging_wall = float(data[12][len("hanging_wall = "):-2])
	return ovb_thickness, salt_thickness, hanging_wall
```

### Block 3: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def run(formulation):
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_overburden_coarse")
	grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 4: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder = os.path.join("output", "case_0", f"{formulation}")
```

### Block 5: Extract region indices

The region-index dictionaries come from `GridHandlerGMSH`. They let the script assign arrays element by element, so salt, overburden, or separate cube regions can share one mesh while keeping distinct properties.

```Python
	# Extract region indices
	ind_salt = grid.region_indices["Salt"]
	ind_ovb = grid.region_indices["Overburden"]
```

### Block 6: Define momentum equation

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

### Block 7: Define solver

The mechanical linear system is handed to PETSc. GMRES is used for the momentum formulations because mixed forms, stabilization, and history-dependent tangents can make the assembled system less friendly to simple symmetric solvers.

```Python
	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("gmres")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)
```

### Block 8: Define material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python
	# Define material properties
	mat = sf.Material(mom_eq.n_elems)
```

### Block 9: Set material density

Density is assigned as an element-wise tensor rather than a scalar so each geological region can carry its own body-force contribution. The same vector is also available to material routines that need density-dependent quantities.

```Python
	# Set material density
	gas_density = 0.082
	salt_density = 2200
	ovb_density = 2800
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

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Create Kelvin-Voigt viscoelastic element
	# eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)
```

### Block 12: Continuation of the script

The Kelvin-Voigt `Viscoelastic` element adds delayed, rate-dependent strain. It is included when the example needs short-term transient deformation beyond purely elastic response.

```Python
	eta = to.zeros(mom_eq.n_elems)
	eta[ind_salt] = 105e11
	eta[ind_ovb] = 105e21
	kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")
```

### Block 13: Create dislocation creep

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
	# Create dislocation creep
	A = to.zeros(mom_eq.n_elems)
	A[ind_salt] = 1.9e-20
	A[ind_ovb] = 0.0
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_ds = sf.DislocationCreep(A, Q, n, "ds_creep")
```

### Block 14: Create pressure solution creep

`PressureSolutionCreep` represents another salt-creep pathway, controlled by stress, temperature, activation energy, and grain-size parameter `d`. It is combined with dislocation creep through the `Material` non-elastic element list.

```Python
	# Create pressure solution creep
	A = to.zeros(mom_eq.n_elems)
	A[ind_salt] = 1.29e-19
	A[ind_ovb] = 0.0
	Q = 13184*to.ones(mom_eq.n_elems)
	d = 0.01*to.ones(mom_eq.n_elems)
	creep_ps = sf.PressureSolutionCreep(A, d, Q, "ps_creep")
```

### Block 15: Thermo-elastic element

The thermoelastic element converts temperature changes into thermal strain, `alpha * deltaT * I`. Adding it to the material is what lets heat diffusion produce mechanical stress and deformation in coupled examples.

```Python
	# Thermo-elastic element
	alpha = to.zeros(mom_eq.n_elems)
	alpha[ind_salt] = 44e-6
	alpha[ind_ovb] = 0.0
	thermo = sf.Thermoelastic(alpha, "thermo")
```

### Block 16: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms, thermal strain coupling with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_thermoelastic(thermo)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_ds)
	mat.add_to_non_elastic(creep_ps)
```

### Block 17: Set constitutive model

Attaching the material to the momentum equation initializes stiffness-related fields and makes the constitutive mechanisms available to the mechanical solver.

```Python
	# Set constitutive model
	mom_eq.set_material(mat)
```

### Block 18: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)
```

### Block 19: Set initial temperature field

The temperature function is evaluated at element centroids. This is the right representation for mechanical constitutive laws because creep and thermoelastic strain are stored per element.

```Python
	# Set initial temperature field
	km = 1000
	dTdZ = 27/km
	T_top = 273 + 20
	T_field_fun = lambda x,y,z: T_top + dTdZ*(660 - z)
	T0_field = create_field_elems(grid, T_field_fun)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)
```

### Block 20: Time settings for equilibrium stage

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
	# Time settings for equilibrium stage
	tc_eq = sf.TimeControllerParabolic(n_time_steps=20, initial_time=0.0, final_time=10, time_unit="day")
	# tc_eq = sf.TimeController(dt=0.1, final_time=5, initial_time=0.0, time_unit="day")
```

### Block 21: Boundary conditions

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Boundary conditions
	time_values = [0*hour,  1*hour]
	nt = len(time_values)
```

### Block 22: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_west_salt = momBC.DirichletBC(boundary_name="West_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_west_ovb = momBC.DirichletBC(boundary_name = "West_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
```

### Block 23: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_east_salt = momBC.DirichletBC(boundary_name="East_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_east_ovb = momBC.DirichletBC(boundary_name = "East_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
```

### Block 24: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_bottom = momBC.DirichletBC(boundary_name="Bottom", component=2, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
```

### Block 25: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_south_salt = momBC.DirichletBC(boundary_name="South_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_south_ovb = momBC.DirichletBC(boundary_name="South_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
```

### Block 26: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_north_salt = momBC.DirichletBC(boundary_name="North_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
	bc_north_ovb = momBC.DirichletBC(boundary_name="North_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
```

### Block 27: Extract geometry dimensions

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Extract geometry dimensions
	Lx = grid.Lx
	Ly = grid.Ly
	Lz = grid.Lz
	z_surface = 0.0
```

### Block 28: Continuation of the script

The mesh is loaded through `GridHandlerGMSH`, which does more than read coordinates and cells. It also exposes named facet tags for boundary conditions, region tags for material assignment, mesh dimensions, cell volumes, and smoothing operators used when element stresses are reported at nodes or smoothed back to elements.

```Python
	g = 9.81
	ovb_thickness, salt_thickness, hanging_wall = get_geometry_parameters(grid_path)
	cavern_roof = ovb_thickness + hanging_wall
	p_roof = 0 + salt_density*g*hanging_wall + ovb_density*g*ovb_thickness
```

### Block 29: Pressure at the top of the salt layer (bottom of overburden)

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Pressure at the top of the salt layer (bottom of overburden)
	p_top = ovb_density*g*ovb_thickness
```

### Block 30: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = z_surface,
						values = [0*MPa, 0*MPa],
						time_values = [0*day,  10*day],
						g = g_vec[2])
```

### Block 31: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = [0.8*p_roof, 0.8*p_roof],
						time_values = [0*day,  tc_eq.t_final],
						g = g_vec[2])
```

### Block 32: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_equilibrium = momBC.BcHandler(mom_eq)
	bc_equilibrium.add_boundary_condition(bc_west_salt)
	bc_equilibrium.add_boundary_condition(bc_west_ovb)
	bc_equilibrium.add_boundary_condition(bc_east_salt)
	bc_equilibrium.add_boundary_condition(bc_east_ovb)
	bc_equilibrium.add_boundary_condition(bc_bottom)
	bc_equilibrium.add_boundary_condition(bc_south_salt)
	bc_equilibrium.add_boundary_condition(bc_south_ovb)
	bc_equilibrium.add_boundary_condition(bc_north_salt)
	bc_equilibrium.add_boundary_condition(bc_north_ovb)
	bc_equilibrium.add_boundary_condition(bc_top)
	bc_equilibrium.add_boundary_condition(bc_cavern)
```

### Block 33: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_equilibrium)
```

### Block 34: Equilibrium output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Equilibrium output folder
	ouput_folder_equilibrium = os.path.join(output_folder, "equilibrium")
```

### Block 35: Print output folder

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(ouput_folder_equilibrium)
```

### Block 36: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(ouput_folder_equilibrium)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]
```

### Block 37: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
	# Define simulator
	sim = sf.Simulator_M(mom_eq, tc_eq, outputs, True)
	sim.run()
```

### Block 38: Time settings for operation stage

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python





	# Time settings for operation stage
	tc_op = sf.TimeController(dt=0.5, initial_time=0.0, final_time=240, time_unit="day")
```

### Block 39: Define heat diffusion equation

`HeatDiffusion` allocates the scalar temperature space, cell-wise material fields, UFL measures, and temperature storage used by the thermal solve. The example creates it before assigning material properties because `set_material` copies those properties into DG0 finite-element fields.

```Python
	# Define heat diffusion equation
	heat_eq = sf.HeatDiffusion(grid)
```

### Block 40: Define solver

The thermal linear system is solved with PETSc. Conjugate gradients match the symmetric diffusion operator, and the additive Schwarz preconditioner keeps the setup usable in both serial and parallel runs.

```Python
	# Define solver
	solver_heat = PETSc.KSP().create(grid.mesh.comm)
	solver_heat.setType("cg")
	solver_heat.getPC().setType("asm")
	solver_heat.setTolerances(rtol=1e-12, max_it=100)
	heat_eq.set_solver(solver_heat)
```

### Block 41: Set specific heat capacity

Specific heat capacity controls thermal inertia. SafeInCave stores it on the material and copies it into a DG0 field inside `HeatDiffusion`, where it multiplies the transient temperature term.

```Python
	# Set specific heat capacity
	cp = 850*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_specific_heat_capacity(cp)
```

### Block 42: Set thermal conductivity

Thermal conductivity sets the strength of the diffusion operator. Defining it per element keeps the example ready for heterogeneous thermal properties even when the current values are uniform.

```Python
	# Set thermal conductivity
	k = 7*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_thermal_conductivity(k)
```

### Block 43: Set material properties to heat_equation

Attaching the material to `HeatDiffusion` copies conductivity, density, and heat capacity into finite-element fields used by the thermal variational form.

```Python
	# Set material properties to heat_equation
	heat_eq.set_material(mat)
```

### Block 44: Set initial temperature

The initial temperature initializes both current and previous thermal states inside `HeatDiffusion`. Without this step the transient solve would not have a physically meaningful starting field.

```Python
	# Set initial temperature
	T0_field_nodes = create_field_nodes(grid, T_field_fun)
	heat_eq.set_initial_T(T0_field_nodes)
```

### Block 45: Define boundary conditions for heat diffusion

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Define boundary conditions for heat diffusion
	time_values = [tc_op.t_initial, tc_op.t_final]
	nt = len(time_values)
```

### Block 46: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler = heatBC.BcHandler(heat_eq)
```

### Block 47: Continuation of the script

This Dirichlet condition fixes temperature on a named boundary. The values are time-interpolated by the heat `BcHandler`, which turns the user schedule into DOLFINx essential boundary conditions at each step.

```Python
	bc_top = heatBC.DirichletBC("Top", nt*[T_top], time_values)
	bc_bottom = heatBC.NeumannBC("Bottom", nt*[dTdZ], time_values)
	bc_east_salt = heatBC.NeumannBC("East_salt", nt*[0.0], time_values)
	bc_east_ovb = heatBC.NeumannBC("East_ovb", nt*[0.0], time_values)
	bc_west_salt = heatBC.NeumannBC("West_salt", nt*[0.0], time_values)
	bc_west_ovb = heatBC.NeumannBC("West_ovb", nt*[0.0], time_values)
	bc_south_salt = heatBC.NeumannBC("South_salt", nt*[0.0], time_values)
	bc_south_ovb = heatBC.NeumannBC("South_ovb", nt*[0.0], time_values)
	bc_north_salt = heatBC.NeumannBC("North_salt", nt*[0.0], time_values)
	bc_north_ovb = heatBC.NeumannBC("North_ovb", nt*[0.0], time_values)
```

### Block 48: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler.add_boundary_condition(bc_top)
	bc_handler.add_boundary_condition(bc_bottom)
	bc_handler.add_boundary_condition(bc_east_salt)
	bc_handler.add_boundary_condition(bc_east_ovb)
	bc_handler.add_boundary_condition(bc_west_salt)
	bc_handler.add_boundary_condition(bc_west_ovb)
	bc_handler.add_boundary_condition(bc_south_salt)
	bc_handler.add_boundary_condition(bc_south_ovb)
	bc_handler.add_boundary_condition(bc_north_salt)
	bc_handler.add_boundary_condition(bc_north_ovb)
```

### Block 49: Continuation of the script

The Robin condition represents convective heat transfer between the rock boundary and a fluid temperature history. SafeInCave contributes both the boundary stiffness term and the ambient-temperature load term during heat assembly.

```Python
	T_gas = T_top
	h_conv = 5.0
	bc_cavern = heatBC.RobinBC("Cavern", nt*[T_gas], h_conv, time_values)
	bc_handler.add_boundary_condition(bc_cavern)
```

### Block 50: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	heat_eq.set_boundary_conditions(bc_handler)
```

### Block 51: Set operation stage settings for momentum equation

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python




	# Set operation stage settings for momentum equation
```

### Block 52: Boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	# Boundary conditions
	bc_west_salt = momBC.DirichletBC(boundary_name="West_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_west_ovb = momBC.DirichletBC(boundary_name = "West_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
```

### Block 53: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_east_salt = momBC.DirichletBC(boundary_name="East_salt", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_east_ovb = momBC.DirichletBC(boundary_name = "East_ovb", component=0, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
```

### Block 54: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_bottom = momBC.DirichletBC(boundary_name="Bottom", component=2, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
```

### Block 55: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_south_salt = momBC.DirichletBC(boundary_name="South_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_south_ovb = momBC.DirichletBC(boundary_name="South_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
```

### Block 56: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_north_salt = momBC.DirichletBC(boundary_name="North_salt", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
	bc_north_ovb = momBC.DirichletBC(boundary_name="North_ovb", component=1, values=[0.0, 0.0], time_values=[0.0, tc_op.t_final])
```

### Block 57: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_top = momBC.NeumannBC(boundary_name="Top", direction=2, density=0.0, ref_pos=z_surface, values=[0, 0], time_values=[0, tc_op.t_final], g=g_vec[2])
```

### Block 58: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	p_values = 3*[0.8*p_roof, 0.8*p_roof, 0.2*p_roof, 0.2*p_roof] + [0.8*p_roof]
	t_values = [20*day*i for i in range(13)]
```

### Block 59: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = cavern_roof,
						values = p_values,
						time_values = t_values,
						g = g_vec[2])
```

### Block 60: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python

	bc_operation = momBC.BcHandler(mom_eq)
	bc_operation.add_boundary_condition(bc_west_salt)
	bc_operation.add_boundary_condition(bc_west_ovb)
	bc_operation.add_boundary_condition(bc_east_salt)
	bc_operation.add_boundary_condition(bc_east_ovb)
	bc_operation.add_boundary_condition(bc_bottom)
	bc_operation.add_boundary_condition(bc_south_salt)
	bc_operation.add_boundary_condition(bc_south_ovb)
	bc_operation.add_boundary_condition(bc_north_salt)
	bc_operation.add_boundary_condition(bc_north_ovb)
	bc_operation.add_boundary_condition(bc_top)
	bc_operation.add_boundary_condition(bc_cavern)
```

### Block 61: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_operation)
```

### Block 62: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder_operation = os.path.join(output_folder, "operation")
```

### Block 63: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder_operation)
		sys.stdout.flush()
```

### Block 64: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder_operation)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
```

### Block 65: Continuation of the script

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder_operation)
	output_heat.add_output_field("T", "Temperature (K)")
```

### Block 66: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	outputs = [output_mom, output_heat]
```

### Block 67: Define simulator

`Simulator_TM` couples heat diffusion and mechanics without solving a cavern thermodynamic model. After each heat step, element temperatures are transferred to the mechanical constitutive update.

```Python
	# Define simulator
	sim = sf.Simulator_TM(mom_eq, heat_eq, tc_op, outputs, False)
	sim.run()
```

### Block 68: Function `main`

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
	run("P1")
	run("P1P1")
	run("P1P1_Stab")
```

### Block 69: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python

if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: reads the cavern-overburden mesh.
- `safeincave.LinearMomentum` and `safeincave.LinearMomentumMixed`: define mechanical formulations for the cavern.
- `safeincave.HeatDiffusion`: defines the thermal problem for the operation stage.
- `safeincave.Material`: assigns salt and overburden properties by region.
- `safeincave.Spring`, `safeincave.Viscoelastic`, `safeincave.DislocationCreep`, `safeincave.PressureSolutionCreep`, and `safeincave.Thermoelastic`: define the constitutive behavior.
- `safeincave.Utils.create_field_elems` and `safeincave.Utils.create_field_nodes`: create geothermal temperature fields.
- `safeincave.Simulator_M` and `safeincave.Simulator_TM`: run the mechanical equilibrium stage and coupled operation stage.

## Running the Example

From the repository root, run:

```bash
cd examples/thermomechanics/2_cavern
python main.py
```

The script writes equilibrium and operation outputs under `output/case_0`. After completion, run:

```bash
python plot_results.py
```

## Conclusion

This example shows how to combine an initial geostatic calculation with a later coupled thermomechanical operation. It is a template for cavern storage analyses that need separate initialization and operation phases.
