# Example 1: Cube Heat Diffusion

## Introduction

This example solves transient heat diffusion in a cube. It applies a fixed temperature on one face and convective heat exchange on the opposite face, starting from a uniform initial temperature field.

## Problem description

This section walks through `examples/thermal/1_cube/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: thermal boundary conditions, PETSc linear solvers. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
import safeincave.HeatBC as heatBC
from petsc4py import PETSc
from mpi4py import MPI
import torch as to
import os
import sys
```

### Block 2: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cube")
	grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 3: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder = os.path.join("output", "case_0")
```

### Block 4: Time settings for equilibrium stage

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
	# Time settings for equilibrium stage
	t_control = sf.TimeControllerParabolic(n_time_steps=50, initial_time=0.0, final_time=5, time_unit="day")
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
							values = [273, 273],
							time_values = [t_control.t_initial, t_control.t_final])
```

### Block 14: Continuation of the script

The Robin condition represents convective heat transfer between the rock boundary and a fluid temperature history. SafeInCave contributes both the boundary stiffness term and the ambient-temperature load term during heat assembly.

```Python
	bc_west = heatBC.RobinBC(boundary_name = "WEST",
							values = [273, 273],
							h = 5.0,
							time_values = [t_control.t_initial, t_control.t_final])
```

### Block 15: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler = heatBC.BcHandler(heat_eq)
	bc_handler.add_boundary_condition(bc_east)
	bc_handler.add_boundary_condition(bc_west)
```

### Block 16: Add boundary condition to heat equation

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Add boundary condition to heat equation
	heat_eq.set_boundary_conditions(bc_handler)
```

### Block 17: Set initial temperature field

The initial temperature initializes both current and previous thermal states inside `HeatDiffusion`. Without this step the transient solve would not have a physically meaningful starting field.

```Python
	# Set initial temperature field
	T0_field = 293*to.ones(heat_eq.n_elems, dtype=to.float64)
	heat_eq.set_initial_T(T0_field)
```

### Block 18: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder)
	output_heat.add_output_field("T", "Temperature (K)")
	outputs = [output_heat]
```

### Block 19: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder)
		sys.stdout.flush()
```

### Block 20: Define simulator

`Simulator_T` advances the heat equation and writes the thermal fields. It is the appropriate driver when the example is only concerned with temperature evolution.

```Python
	# Define simulator
	sim = sf.Simulator_T(heat_eq, t_control, outputs, True)
	sim.run()
```

### Block 21: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python


if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: loads the cube mesh from `grids/cube`.
- `safeincave.HeatDiffusion`: defines the heat diffusion equation.
- `safeincave.Material`: stores density, specific heat capacity, and thermal conductivity.
- `safeincave.HeatBC.DirichletBC`: prescribes a fixed temperature on the `EAST` boundary.
- `safeincave.HeatBC.RobinBC`: applies convective heat exchange on the `WEST` boundary.
- `safeincave.TimeControllerParabolic`: creates the parabolic time discretization.
- `safeincave.SaveFields`: saves the temperature field.
- `safeincave.Simulator_T`: runs the thermal simulation.

## Running the Example

From the repository root, run:

```bash
cd examples/thermal/1_cube
python main.py
```

The simulation writes temperature results to `output/case_0`. To inspect the results with the provided plotting script, run:

```bash
python plot_results.py
```

## Conclusion

This example demonstrates the minimal workflow for a thermal SafeInCave simulation: load a mesh, define thermal material properties, assign thermal boundary conditions, save temperature, and run `Simulator_T`.
