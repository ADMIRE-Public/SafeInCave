# Example 1: Thermoelastic Cube

## Introduction

This example runs a coupled thermomechanical simulation on a cube with two material regions. A cooling boundary condition produces a temperature field, and the mechanical solver computes the resulting thermoelastic stresses for three formulations: `P1`, `P1P1`, and stabilized `P1P1`.

## Problem description

This section walks through `examples/thermomechanics/1_cube/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, thermal boundary conditions, PETSc linear solvers. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
import safeincave.Utils as ut
import safeincave.HeatBC as heatBC
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import os
import sys
```

### Block 2: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python


def run(formulation):
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube_regions")
	grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 3: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder = os.path.join("output", "case_0", f"{formulation}")
```

### Block 4: Time settings for equilibrium stage

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
	# Time settings for equilibrium stage
	t_control = sf.TimeControllerParabolic(n_time_steps=100, initial_time=0.0, final_time=10, time_unit="day")
```

### Block 5: Define equation

`HeatDiffusion` allocates the scalar temperature space, cell-wise material fields, UFL measures, and temperature storage used by the thermal solve. The example creates it before assigning material properties because `set_material` copies those properties into DG0 finite-element fields.

```Python
	# Define equation
	heat_eq = sf.HeatDiffusion(grid)
```

### Block 6: Define solver

The thermal linear system is solved with PETSc. Conjugate gradients match the symmetric diffusion operator, and the additive Schwarz preconditioner keeps the setup usable in both serial and parallel runs.

```Python
	# Define solver
	solver_heat = PETSc.KSP().create(grid.mesh.comm)
	solver_heat.setType("cg")
	solver_heat.getPC().setType("asm")
	solver_heat.setTolerances(rtol=1e-12, max_it=100)
	heat_eq.set_solver(solver_heat)
```

### Block 7: Build material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python
	# Build material properties
	mat = sf.Material(heat_eq.n_elems)
```

### Block 8: Set material density

The density vector is stored in the `Material` object. Mechanics uses it when building body forces, while thermal examples keep it together with heat capacity for the transient heat-storage term.

```Python
	# Set material density
	rho = 2000.0*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)
```

### Block 9: Set specific heat capacity

Specific heat capacity controls thermal inertia. SafeInCave stores it on the material and copies it into a DG0 field inside `HeatDiffusion`, where it multiplies the transient temperature term.

```Python
	# Set specific heat capacity
	cp = 850*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_specific_heat_capacity(cp)
```

### Block 10: Set thermal conductivity

Thermal conductivity sets the strength of the diffusion operator. Defining it per element keeps the example ready for heterogeneous thermal properties even when the current values are uniform.

```Python
	# Set thermal conductivity
	k = 7*to.ones(heat_eq.n_elems, dtype=to.float64)
	mat.set_thermal_conductivity(k)
```

### Block 11: Set material properties to heat_equation

Attaching the material to `HeatDiffusion` copies conductivity, density, and heat capacity into finite-element fields used by the thermal variational form.

```Python
	# Set material properties to heat_equation
	heat_eq.set_material(mat)
```

### Block 12: Define boundary conditions for heat diffusion

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Define boundary conditions for heat diffusion
	time_values = [t_control.t_initial, t_control.t_final]
	nt = len(time_values)
```

### Block 13: Continuation of the script

This Dirichlet condition fixes temperature on a named boundary. The values are time-interpolated by the heat `BcHandler`, which turns the user schedule into DOLFINx essential boundary conditions at each step.

```Python
	bc_east = heatBC.DirichletBC(boundary_name = "EAST",
							values = nt*[273 + 1.0],
							time_values = time_values)
```

### Block 14: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler = heatBC.BcHandler(heat_eq)
	bc_handler.add_boundary_condition(bc_east)
	heat_eq.set_boundary_conditions(bc_handler)
```

### Block 15: Set initial temperature field

The initial temperature initializes both current and previous thermal states inside `HeatDiffusion`. Without this step the transient solve would not have a physically meaningful starting field.

```Python
	# Set initial temperature field
	fun = lambda x, y, z: 273 + 20
	T0_field = ut.create_field_nodes(heat_eq.grid, fun)
	heat_eq.set_initial_T(T0_field)
```

### Block 16: Define momentum equation

The mixed momentum formulation is selected when pressure-related stability or volumetric behavior matters. The `stab_scaling` parameter controls the pressure-stabilization contribution, while `theta` sets the time-integration convention used by inelastic updates.

```Python


	# Define momentum equation
	theta = 0.5
	if formulation == "P1":
		mom_eq = sf.LinearMomentum(grid, theta=theta)
	elif formulation == "P1P1":
		mom_eq = sf.LinearMomentumMixed(grid, theta=theta, stab_scaling=0.0)
	elif formulation == "P1P1_Stab":
		mom_eq = sf.LinearMomentumMixed(grid, theta=theta, stab_scaling=1.0)
```

### Block 17: Define solver

The mechanical linear system is handed to PETSc. GMRES is used for the momentum formulations because mixed forms, stabilization, and history-dependent tangents can make the assembled system less friendly to simple symmetric solvers.

```Python
	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("gmres")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)
```

### Block 18: Constitutive model

The `Spring` element supplies the isotropic elastic stiffness matrix used as the backbone of the mechanical model. Other inelastic mechanisms evolve around this elastic response.

```Python
	# Constitutive model
	E = 102e9*to.ones(mom_eq.n_elems)
	nu = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E, nu, "spring")
```

### Block 19: Extract region indices

The region-index dictionaries come from `GridHandlerGMSH`. They let the script assign arrays element by element, so salt, overburden, or separate cube regions can share one mesh while keeping distinct properties.

```Python
	# Extract region indices
	omega_A = grid.region_indices["OMEGA_A"]
	omega_B = grid.region_indices["OMEGA_B"]
```

### Block 20: Thermo-elastic element

The thermoelastic element converts temperature changes into thermal strain, `alpha * deltaT * I`. Adding it to the material is what lets heat diffusion produce mechanical stress and deformation in coupled examples.

```Python
	# Thermo-elastic element
	alpha = to.zeros(mom_eq.n_elems)
	alpha[omega_A] = 44e-6
	alpha[omega_B] = 74e-6
	# alpha = 44e-6*to.ones(mom_eq.n_elems)
	thermo = sf.Thermoelastic(alpha, "thermo")
```

### Block 21: Create constitutive model

These calls register elastic stiffness, thermal strain coupling with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_thermoelastic(thermo)
```

### Block 22: Set constitutive model

Attaching the material to the momentum equation initializes stiffness-related fields and makes the constitutive mechanisms available to the mechanical solver.

```Python
	# Set constitutive model
	mom_eq.set_material(mat)
```

### Block 23: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
	# Set body forces
	g = -9.81
	g_vec = [0.0, 0.0, g]
	mom_eq.build_body_force(g_vec)
```

### Block 24: Boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	# Boundary conditions
	bc_west_2 = momBC.DirichletBC(boundary_name = "WEST",
							component = 2,
							values = nt*[0.0],
							time_values = time_values)
```

### Block 25: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_west_1 = momBC.DirichletBC(boundary_name = "WEST",
							component = 1,
							values = nt*[0.0],
							time_values = time_values)
```

### Block 26: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_west_0 = momBC.DirichletBC(boundary_name = "WEST",
						  component = 0,
						  values = nt*[0.0],
						  time_values = time_values)
```

### Block 27: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_bottom = momBC.DirichletBC(boundary_name = "BOTTOM",
						  component = 2,
						  values = nt*[0.0],
						  time_values = time_values)
```

### Block 28: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler = momBC.BcHandler(mom_eq)
	bc_handler.add_boundary_condition(bc_west_0)
	bc_handler.add_boundary_condition(bc_west_1)
	bc_handler.add_boundary_condition(bc_west_2)
	bc_handler.add_boundary_condition(bc_bottom)
```

### Block 29: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_handler)
```

### Block 30: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("sig", "Stress (Pa)")
	output_mom.add_output_field("p_nodes", "Mean stress (Pa)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_nodes", "Von Mises stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
```

### Block 31: Continuation of the script

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder)
	output_heat.add_output_field("T", "Temperature (K)")
```

### Block 32: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	outputs = [output_mom, output_heat]
```

### Block 33: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)
		sys.stdout.flush()
```

### Block 34: Define simulator

`Simulator_TM` couples heat diffusion and mechanics without solving a cavern thermodynamic model. After each heat step, element temperatures are transferred to the mechanical constitutive update.

```Python
	# Define simulator
	sim = sf.Simulator_TM(mom_eq, heat_eq, t_control, outputs, True)
	sim.run()
```

### Block 35: Function `main`

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
	run("P1")
	run("P1P1")
	run("P1P1_Stab")
```

### Block 36: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python


if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: loads the `cube_regions` mesh and region labels.
- `safeincave.HeatDiffusion`: solves transient heat diffusion.
- `safeincave.LinearMomentum` and `safeincave.LinearMomentumMixed`: define the mechanical formulations.
- `safeincave.Material`: stores thermal and mechanical material properties.
- `safeincave.Spring` and `safeincave.Thermoelastic`: define elastic stiffness and thermal expansion.
- `safeincave.HeatBC.DirichletBC`: applies a fixed temperature boundary.
- `safeincave.MomentumBC.DirichletBC`: constrains selected displacement components.
- `safeincave.Simulator_TM`: runs the coupled thermal-mechanical problem.

## Running the Example

From the repository root, run:

```bash
cd examples/thermomechanics/1_cube
python main.py
```

The script runs all three mechanical formulations and writes separate result folders under `output/case_0`. To inspect the fields, run:

```bash
python plot_results.py
```

## Conclusion

This example demonstrates the basic coupled workflow in SafeInCave and shows how thermal expansion coefficients can vary by mesh region. It also provides a compact comparison setup for the available mechanical formulations.
