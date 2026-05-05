# Cube with 2 regions

## Goal:

1. Define material properties in different regions

## Problem description

This problem consists of a cube divided in two regions: OMEGA_A and OMEGA_B, as illustrated in [](#fig-cube-geom). The cube is subjected to a constant confining pressure of 5 MPa on faces EAST and NORTH, and a constant axial stress of 8 MPa on face TOP. The remaining faces are prevented from normal displacements.

![(a) Boundary names; (b) Region names; (c) Axial load and confining pressure.](../../images/2_cube_regions_geom.png){#fig-cube-geom width="100%"}

As illustrated in [](#fig-cube-model), the constitutive model consists of an elastic element (spring) and a viscoelastic (kelvin) element. Different values for the material properties of these two elements are assigned to the two regions of the domain, as indicated in the table of [](#fig-cube-model).

![Viscoelastic constitutive modeland material parameters for regions OMEGA_A and OMEGA_B.](../../images/2_cube_regions_model.png){#fig-cube-model width="40%"}

### Complete annotated `main.py`

The listing below follows `examples/mechanics/2_cube_regions/main.py` in order and includes all non-empty source lines. The original overview above describes the physical problem; this listing documents the complete executable workflow and explains why each SafeInCave object is introduced.

#### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, PETSc linear solvers. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
from safeincave import *
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from mpi4py import MPI
from petsc4py import PETSc
import torch as to
import os
```

#### Block 2: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def run(formulation):
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube_regions")
	grid = GridHandlerGMSH("geom", grid_path)
```

#### Block 3: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder = os.path.join("output", "case_0", f"{formulation}")
```

#### Block 4: Time settings for equilibrium stage

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python
	# Time settings for equilibrium stage
	unit = "hour"
	t_0 = 0.0
	dt = 0.01
	t_final = 1
	t_control = TimeController(dt=dt, initial_time=t_0, final_time=t_final, time_unit=unit)
```

#### Block 5: Define momentum equation

The mixed momentum formulation is selected when pressure-related stability or volumetric behavior matters. The `stab_scaling` parameter controls the pressure-stabilization contribution, while `theta` sets the time-integration convention used by inelastic updates.

```Python
	# Define momentum equation
	theta = 0.5
	if formulation == "P1":
		mom_eq = LinearMomentum(grid, theta=theta)
	elif formulation == "P1P1":
		mom_eq = LinearMomentumMixed(grid, theta=theta, stab_scaling=0.0)
	elif formulation == "P1P1_Stab":
		mom_eq = LinearMomentumMixed(grid, theta=theta, stab_scaling=1.0)
```

#### Block 6: Define solver

The mechanical linear system is handed to PETSc. GMRES is used for the momentum formulations because mixed forms, stabilization, and history-dependent tangents can make the assembled system less friendly to simple symmetric solvers.

```Python

	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("gmres")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)
```

#### Block 7: Define material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python
	# Define material properties
	mat = Material(mom_eq.n_elems)
```

#### Block 8: Set material density

The density vector is stored in the `Material` object. Mechanics uses it when building body forces, while thermal examples keep it together with heat capacity for the transient heat-storage term.

```Python
	# Set material density
	rho = 0.0*to.ones(mom_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)
```

#### Block 9: Extract region indices

The region-index dictionaries come from `GridHandlerGMSH`. They let the script assign arrays element by element, so salt, overburden, or separate cube regions can share one mesh while keeping distinct properties.

```Python
	# Extract region indices
	omega_A = grid.region_indices["OMEGA_A"]
	omega_B = grid.region_indices["OMEGA_B"]
```

#### Block 10: Constitutive model

The spring is the elastic part of the constitutive model. Here its Young modulus and Poisson ratio are assembled as element-wise arrays, which is why different mesh regions can receive different elastic stiffness.

```Python
	# Constitutive model
	E0 = to.zeros(mom_eq.n_elems)
	nu0 = to.zeros(mom_eq.n_elems)
	E0[omega_A] = 8*ut.GPa
	E0[omega_B] = 10*ut.GPa
	nu0[omega_A] = 0.2
	nu0[omega_B] = 0.3
	spring_0 = Spring(E0, nu0, "spring")
```

#### Block 11: Create Kelvin-Voigt viscoelastic element

The Kelvin-Voigt `Viscoelastic` element adds delayed, rate-dependent strain. It is included when the example needs short-term transient deformation beyond purely elastic response.

```Python
	# Create Kelvin-Voigt viscoelastic element
	eta = to.zeros(mom_eq.n_elems)
	E1 = to.zeros(mom_eq.n_elems)
	nu1 = to.zeros(mom_eq.n_elems)
	eta[omega_A] = 105e11
	eta[omega_B] = 38e11
	E1[omega_A] = 8*ut.GPa
	E1[omega_B] = 5*ut.GPa
	nu1[omega_A] = 0.35
	nu1[omega_B] = 0.28
	kelvin = Viscoelastic(eta, E1, nu1, "kelvin")
```

#### Block 12: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)
```

#### Block 13: Set constitutive model

Attaching the material to the momentum equation initializes stiffness-related fields and makes the constitutive mechanisms available to the mechanical solver.

```Python
	# Set constitutive model
	mom_eq.set_material(mat)
```

#### Block 14: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)
```

#### Block 15: Set initial temperature field

The temperature arrays define the reference and current thermal state used by temperature-dependent mechanical mechanisms.

```Python
	# Set initial temperature field
	T0_field = 298*to.ones(mom_eq.n_elems)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)
```

#### Block 16: Boundary conditions

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Boundary conditions
	time_values = [0*ut.hour,  1*ut.hour]
	nt = len(time_values)
```

#### Block 17: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_west = momBC.DirichletBC(boundary_name = "WEST",
							component = 0,
							values = [0.0, 0.0],
							time_values = [0.0, t_control.t_final])
```

#### Block 18: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_bottom = momBC.DirichletBC(boundary_name = "BOTTOM",
						  component = 2,
						  values = [0.0, 0.0],
						  time_values = [0.0, t_control.t_final])
```

#### Block 19: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_south = momBC.DirichletBC(boundary_name = "SOUTH",
						  component = 1,
						  values = [0.0, 0.0],
						  time_values = [0.0, t_control.t_final])
```

#### Block 20: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_east = momBC.NeumannBC(boundary_name = "EAST",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [5.0*ut.MPa, 5.0*ut.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])
```

#### Block 21: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_north = momBC.NeumannBC(boundary_name = "NORTH",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [5.0*ut.MPa, 5.0*ut.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])
```

#### Block 22: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_top = momBC.NeumannBC(boundary_name = "TOP",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [8.0*ut.MPa, 8.0*ut.MPa],
						time_values = [0.0,           t_control.t_final],
						g = g_vec[2])
```

#### Block 23: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler = momBC.BcHandler(mom_eq)
	bc_handler.add_boundary_condition(bc_west)
	bc_handler.add_boundary_condition(bc_bottom)
	bc_handler.add_boundary_condition(bc_south)
	bc_handler.add_boundary_condition(bc_east)
	bc_handler.add_boundary_condition(bc_north)
	bc_handler.add_boundary_condition(bc_top)
```

#### Block 24: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_handler)
```

#### Block 25: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]
```

#### Block 26: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)
```

#### Block 27: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
	# Define simulator
	sim = Simulator_M(mom_eq, t_control, outputs, compute_elastic_response=True)
	sim.run()
```

#### Block 28: run("P1")

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
	# run("P1")
	# run("P1P1")
	run("P1P1_Stab")
```

#### Block 29: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python


if __name__ == '__main__':
	main()
```

## Implementation

As usual, import relevant packages.

```Python
from safeincave import *
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
import torch as to
import os
import time
```


Define grid path and create grid object.

```Python
grid_path = os.path.join("..", "..", "..", "grids", "cube_regions")
grid = GridHandlerGMSH("geom", grid_path)
```

Define output folder where the simulation results will be saved.

```Python
output_folder = os.path.join("output", "case_0")
```

Create an equally spaced time discretization with time step size of 0.01 hour and final time of 1.0 hour. For this purpose, we use class **TimeController**.

```Python
t_control = TimeController(dt=0.01, initial_time=0.0, final_time=1.0, time_unit="hour")
```

Initialize object for the momentum balance equation (**LinearMomentum**) and choose Crank-Nicolson as a time integration scheme ($\theta=0.5$).

```Python
mom_eq = LinearMomentum(grid, theta=0.5)
```

Define solver for momentum balance equation. Choose Conjugate Gradient as a linear solver with Additive Schwartz preconditioner.

> **Note:**: The reason we use Additive Schwartz Method here is because it works well in series and parallel. For instance, incomplete LU factorization (ILU) works well with serial computations, but not in parallel.

```Python
mom_solver = PETSc.KSP().create(grid.mesh.comm)
mom_solver.setType("bicg")
mom_solver.getPC().setType("asm")
mom_solver.setTolerances(rtol=1e-12, max_it=100)
```

Set solver to momentum equation object.

```Python
mom_eq.set_solver(mom_solver)
```

Initialize **Material** object, which contains all material properties and the constitutive model.

```Python
mat = Material(mom_eq.n_elems)
```

Define and set zero density to eliminate body forces effect.

```Python
rho = 0.0*to.ones(mom_eq.n_elems, dtype=to.float64)
mat.set_density(rho)
```

Extract lists of indices belonging to regions OMEGA_A and OMEGA_B. Notice that attribute *region_indices* is a dictionary with as many keys as the number of regions in the mesh.

```Python
omega_A = grid.region_indices["OMEGA_A"]
omega_B = grid.region_indices["OMEGA_B"]
```

Create spring element. First, create a vector of zeros for property $E_0$, and then assign 8 GPa to all elements in OMEGA_A, and 10 GPa to elements in OMEGA_B. Do the same for Poisson's ratio, $\nu_0$. Finally, pass these arguments to class **Spring**, together with a given name '*spring*'.

```Python
E0 = to.zeros(mom_eq.n_elems)
E0[omega_A] = 8*ut.GPa
E0[omega_B] = 10*ut.GPa
nu0 = to.zeros(mom_eq.n_elems)
nu0[omega_A] = 0.2
nu0[omega_B] = 0.3
spring_0 = Spring(E0, nu0, "spring")
```

Create the viscoelastic (i.e., Kelvin-Voigt) element. Assign different values for $\eta_1$, $E_1$, and $\nu_1$ to regions OMEGA_A and OMEGA_B according to the table in" + f" Fig. {fig_2_cube_region_model}.

```Python
eta = to.zeros(mom_eq.n_elems)
eta[omega_A] = 105e11
eta[omega_B] = 38e11
E1 = to.zeros(mom_eq.n_elems)
E1[omega_A] = 8*ut.GPa
E1[omega_B] = 5*ut.GPa
nu1 = to.zeros(mom_eq.n_elems)
nu1[omega_A] = 0.35
nu1[omega_B] = 0.28
kelvin = Viscoelastic(eta, E1, nu1, "kelvin")
```

Add element *spring* to the elastic list of the material *mat*, and add element *kelvin* to the list of non-elastic elements of *mat*.

```Python
mat.add_to_elastic(spring_0)
mat.add_to_non_elastic(kelvin)
```

Now that the material is completely defined, set it to the momentum equation object.

```Python
mom_eq.set_material(mat)
```

Define gravity acceleration vector and assign it to *mom_eq* so it builds the body force terms.

```Python
g = -9.81
g_vec = [0.0, 0.0, g]
mom_eq.build_body_force(g_vec)
```

Define an uniform temperature field of 298 K and assign it to both initial and current temperature.

>**Note:** For the constitutive model adopted in this example, the temperature field is not particularly important, since neither the spring nor the Kelvin-Voigt elements depend on it. In this manner, it's not really necessary to specify temperature for this case. However, it is a good practice to always specify temperature, as the user might decide later to include, for example, dislocation creep into the constitutive model. If dislocation (or pressure solution) creep is present but temperature is not specified, the program will consider internally temperature to be zero, which will in practice neutralize the creep effect according to Arrhenius law.

```Python
T0_field = 298*to.ones(mom_eq.n_elems)
mom_eq.set_T0(T0_field)
mom_eq.set_T(T0_field)
```

Impose zero normal displacements (Dirichlet boundary conditions) on faces WEST, SOUTH, and BOTTOM, by specifying displacement components 0 (x), 1 (y), and 2 (z), respectively.

```Python
bc_west = momBC.DirichletBC(boundary_name = "WEST",
					component = 0,
					values = [0.0, 0.0],
					time_values = [0.0, t_control.t_final])
bc_south = momBC.DirichletBC(boundary_name = "SOUTH",
					component = 1,
					values = [0.0, 0.0],
					time_values = [0.0, t_control.t_final])
bc_bottom = momBC.DirichletBC(boundary_name = "BOTTOM",
					component = 2,
					values = [0.0, 0.0],
					time_values = [0.0, t_control.t_final])
```

Impose constant and uniform normal loads (Neumann boundary condition) on faces EAST, NORTH, and TOP.

>**Note:** Because the loads are uniform along the boundary, we specify *density* to be zero. Consequently, the arguments *direction* and *ref_pos* are ignored.

```Python
bc_east = momBC.NeumannBC(boundary_name = "EAST",
					direction = 2,
					density = 0.0,
					ref_pos = 0.0,
					values = [5.0*ut.MPa, 5.0*ut.MPa],
					time_values = [0.0, t_control.t_final],
					g = g_vec[2])
bc_north = momBC.NeumannBC(boundary_name = "NORTH",
					direction = 2,
					density = 0.0,
					ref_pos = 0.0,
					values = [5.0*ut.MPa, 5.0*ut.MPa],
					time_values = [0.0, t_control.t_final],
					g = g_vec[2])
bc_top = momBC.NeumannBC(boundary_name = "TOP",
					direction = 2,
					density = 0.0,
					ref_pos = 0.0,
					values = [8.0*ut.MPa, 8.0*ut.MPa],
					time_values = [0.0, t_control.t_final],
					g = g_vec[2])
```

Create a **BcHandler** object and add the above defined boundary conditions to it.

```Python
bc_handler = momBC.BcHandler(mom_eq)
bc_handler.add_boundary_condition(bc_west)
bc_handler.add_boundary_condition(bc_bottom)
bc_handler.add_boundary_condition(bc_south)
bc_handler.add_boundary_condition(bc_east)
bc_handler.add_boundary_condition(bc_north)
bc_handler.add_boundary_condition(bc_top)
```

Set the **BcHandler** object to the momentum balance equation object, *mom_eq*.

```Python
mom_eq.set_boundary_conditions(bc_handler)
```

Choose fields to be saved during the simulation.

> **Note:** As a reminder, the firt argument of *add_output_field* must be a string with the same name of an existing attribute of object *mom_eq*.

```Python
output_mom = SaveFields(mom_eq)
output_mom.set_output_folder(output_folder)
output_mom.add_output_field("u", "Displacement (m)")
output_mom.add_output_field("eps_tot", "Total strain (-)")
output_mom.add_output_field("sig", "Stress (Pa)")
output_mom.add_output_field("p_elems", "Mean stress (Pa)")
output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
outputs = [output_mom]
```

Once objects *mom_eq*, *t_control*, and *outputs* are created, pass them as arguments to the mechanical simulator **Simulator_M**.

```Python
sim = Simulator_M(mom_eq, t_control, outputs, compute_elastic_response=True)
sim.run()
```
