# Example 2: Cavern Heat Diffusion

## Introduction

This example computes heat diffusion around a cavern over a five-year period. The model uses a geothermal gradient as the initial condition, insulated lateral boundaries, and a convective condition at the cavern wall with a time-dependent gas temperature.

## Problem description

This section walks through `examples/thermal/2_cavern/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: thermal boundary conditions, PETSc linear solvers. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
import safeincave.Utils as ut
import safeincave.HeatBC as heatBC
from petsc4py import PETSc
import torch as to
import os
```

### Block 2: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python
def main():
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_regular")
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
	t_control = sf.TimeControllerParabolic(n_time_steps=100, initial_time=0, final_time=5, time_unit="year")
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

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	km = 1000
	dTdZ = 27/km
	T_top = 273 + 20
	h_conv = 5.0
```

### Block 14: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_handler = heatBC.BcHandler(heat_eq)
```

### Block 15: Continuation of the script

This Dirichlet condition fixes temperature on a named boundary. The values are time-interpolated by the heat `BcHandler`, which turns the user schedule into DOLFINx essential boundary conditions at each step.

```Python
	bc_top = heatBC.DirichletBC("Top", nt*[T_top], time_values)
	bc_handler.add_boundary_condition(bc_top)
```

### Block 16: Continuation of the script

This Neumann condition prescribes heat flux. A geothermal gradient at the bottom supplies heat from depth, while zero flux on side boundaries idealizes lateral insulation.

```Python
	bc_bottom = heatBC.NeumannBC("Bottom", nt*[dTdZ], time_values)
	bc_handler.add_boundary_condition(bc_bottom)
```

### Block 17: Continuation of the script

This Neumann condition prescribes heat flux. A geothermal gradient at the bottom supplies heat from depth, while zero flux on side boundaries idealizes lateral insulation.

```Python
	bc_east = heatBC.NeumannBC("East", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_east)
```

### Block 18: Continuation of the script

This Neumann condition prescribes heat flux. A geothermal gradient at the bottom supplies heat from depth, while zero flux on side boundaries idealizes lateral insulation.

```Python
	bc_west = heatBC.NeumannBC("West", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_west)
```

### Block 19: Continuation of the script

This Neumann condition prescribes heat flux. A geothermal gradient at the bottom supplies heat from depth, while zero flux on side boundaries idealizes lateral insulation.

```Python
	bc_south = heatBC.NeumannBC("South", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_south)
```

### Block 20: Continuation of the script

This Neumann condition prescribes heat flux. A geothermal gradient at the bottom supplies heat from depth, while zero flux on side boundaries idealizes lateral insulation.

```Python
	bc_north = heatBC.NeumannBC("North", nt*[0.0], time_values)
	bc_handler.add_boundary_condition(bc_north)
```

### Block 21: Gas temperature schedule over time

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

	# Gas temperature schedule over time
	#      T(K) ^
	#           |
	#       310 |_ _ _  ___________
	#    		|      /|          |
	#    		|    /
	#    		|  /    |          |
	#       303 |/
	#    		|       |          |
	#    		+-------+----------+-----> Time
	#		   t_0    t_f/2       t_f
```

### Block 22: Continuation of the script

The Robin condition represents convective heat transfer between the rock boundary and a fluid temperature history. SafeInCave contributes both the boundary stiffness term and the ambient-temperature load term during heat assembly.

```Python
	bc_cavern = heatBC.RobinBC(
								boundary_name = "Cavern",
								values = [303, 310, 310],
								h = h_conv,
								time_values = [t_control.t_initial, t_control.t_final/2, t_control.t_final]
	)
	bc_handler.add_boundary_condition(bc_cavern)
```

### Block 23: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	heat_eq.set_boundary_conditions(bc_handler)
```

### Block 24: Set initial temperature field

The initial temperature initializes both current and previous thermal states inside `HeatDiffusion`. Without this step the transient solve would not have a physically meaningful starting field.

```Python
	# Set initial temperature field
	fun = lambda x, y, z: T_top - dTdZ*(z - 660)
	T0_field = ut.create_field_nodes(heat_eq.grid, fun)
	heat_eq.set_initial_T(T0_field)
```

### Block 25: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_heat = sf.SaveFields(heat_eq)
	output_heat.set_output_folder(output_folder)
	output_heat.add_output_field("T", "Temperature (K)")
	outputs = [output_heat]
```

### Block 26: Define simulator

`Simulator_T` advances the heat equation and writes the thermal fields. It is the appropriate driver when the example is only concerned with temperature evolution.

```Python
	# Define simulator
	sim = sf.Simulator_T(heat_eq, t_control, outputs, True)
	sim.run()
```

### Block 27: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python

if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: reads the `cavern_regular` mesh and boundary tags.
- `safeincave.HeatDiffusion`: defines the thermal problem.
- `safeincave.Material`: assigns density, heat capacity, and thermal conductivity.
- `safeincave.Utils.create_field_nodes`: builds the initial geothermal temperature field.
- `safeincave.HeatBC.DirichletBC`: fixes the temperature at the top boundary.
- `safeincave.HeatBC.NeumannBC`: applies the basal geothermal gradient and no-flux side boundaries.
- `safeincave.HeatBC.RobinBC`: models heat transfer between the cavern wall and cavern gas.
- `safeincave.Simulator_T`: advances the transient thermal problem.

## Running the Example

From the repository root, run:

```bash
cd examples/thermal/2_cavern
python main.py
```

Results are written to `output/case_0`. The companion script can be used after the run:

```bash
python plot_results.py
```

## Conclusion

This example shows how to impose geological thermal conditions and a cavern-wall heat transfer law. It is a useful starting point for thermal initialization before coupled thermomechanical cavern analyses.
