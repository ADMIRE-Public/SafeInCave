# Example 3: Single Cavern

## Introduction

This example simulates the mechanical response of a salt cavern in two stages. The first stage computes an initial equilibrium state under geostatic loading and cavern pressure. The second stage applies an operating pressure history and extends the material model with a viscoplastic Desai element.

## Problem description

This section walks through `examples/mechanics/3_cavern/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, PETSc linear solvers. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
import safeincave.Utils as ut
import safeincave.MomentumBC as momBC
from petsc4py import PETSc
from mpi4py import MPI
import dolfinx as do
import os
import torch as to
```

### Block 2: Class `LinearMomentumMod`

The subclass adds extra output fields to the standard momentum equation. SafeInCave already computes displacement, strain, stress, pressure, and von Mises stress; this extension exposes selected internal variables from the constitutive model so they can be written by `SaveFields`.

```Python

class LinearMomentumMod(sf.LinearMomentum):
	def __init__(self, grid, theta):
		super().__init__(grid, theta)
```

### Block 3: Function `initialize`

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	def initialize(self) -> None:
		self.C.x.array[:] = to.flatten(self.mat.C)
		self.Fvp = do.fem.Function(self.DG0_1)
		self.alpha = do.fem.Function(self.DG0_1)
		self.eps_vp = do.fem.Function(self.DG0_3x3)
```

### Block 4: Function `run_after_solve`

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python
	def run_after_solve(self):
		try:
			self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[-1].eps_ne_k)
			self.Fvp.x.array[:] = self.mat.elems_ne[-1].Fvp
			self.alpha.x.array[:] = self.mat.elems_ne[-1].alpha
		except:
			pass
```

### Block 5: Class `LinearMomentumMixedMod`

The subclass adds extra output fields to the standard momentum equation. SafeInCave already computes displacement, strain, stress, pressure, and von Mises stress; this extension exposes selected internal variables from the constitutive model so they can be written by `SaveFields`.

```Python

class LinearMomentumMixedMod(sf.LinearMomentumMixed):
	def __init__(self, grid, theta, stab_scaling):
		super().__init__(grid, theta, stab_scaling)
		self.Fvp = do.fem.Function(self.DG0_1)
		self.alpha = do.fem.Function(self.DG0_1)
		self.eps_vp = do.fem.Function(self.DG0_3x3)
```

### Block 6: Function `run_after_solve`

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python
	def run_after_solve(self):
		try:
			self.eps_vp.x.array[:] = to.flatten(self.mat.elems_ne[-1].eps_ne_k)
			self.Fvp.x.array[:] = self.mat.elems_ne[-1].Fvp
			self.alpha.x.array[:] = self.mat.elems_ne[-1].alpha
		except:
			pass
```

### Block 7: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python


def run(formulation):
	# Read grid
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_irregular")
	grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 8: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder = os.path.join("output", "case_0", f"{formulation}")
```

### Block 9: Define momentum equation

The mixed momentum formulation is selected when pressure-related stability or volumetric behavior matters. The `stab_scaling` parameter controls the pressure-stabilization contribution, while `theta` sets the time-integration convention used by inelastic updates.

```Python
	# Define momentum equation
	theta = 0.5
	if formulation == "P1":
		mom_eq = LinearMomentumMod(grid, theta=theta)
	elif formulation == "P1P1":
		mom_eq = LinearMomentumMixedMod(grid, theta=theta, stab_scaling=0.0)
	elif formulation == "P1P1_Stab":
		mom_eq = LinearMomentumMixedMod(grid, theta=theta, stab_scaling=1.0)
```

### Block 10: Define solver

The mechanical linear system is handed to PETSc. GMRES is used for the momentum formulations because mixed forms, stabilization, and history-dependent tangents can make the assembled system less friendly to simple symmetric solvers.

```Python
	# Define solver
	mom_solver = PETSc.KSP().create(grid.mesh.comm)
	mom_solver.setType("gmres")
	mom_solver.getPC().setType("asm")
	mom_solver.setTolerances(rtol=1e-12, max_it=100)
	mom_eq.set_solver(mom_solver)
```

### Block 11: Define material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python
	# Define material properties
	mat = sf.Material(mom_eq.n_elems)
```

### Block 12: Set material density

The density vector is stored in the `Material` object. Mechanics uses it when building body forces, while thermal examples keep it together with heat capacity for the transient heat-storage term.

```Python
	# Set material density
	salt_density = 2000
	rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
	mat.set_density(rho)
```

### Block 13: Constitutive model

The `Spring` element supplies the isotropic elastic stiffness matrix used as the backbone of the mechanical model. Other inelastic mechanisms evolve around this elastic response.

```Python
	# Constitutive model
	E0 = 102*ut.GPa*to.ones(mom_eq.n_elems)
	nu0 = 0.3*to.ones(mom_eq.n_elems)
	spring_0 = sf.Spring(E0, nu0, "spring")
```

### Block 14: Create Kelvin-Voigt viscoelastic element

The Kelvin-Voigt `Viscoelastic` element adds delayed, rate-dependent strain. It is included when the example needs short-term transient deformation beyond purely elastic response.

```Python
	# Create Kelvin-Voigt viscoelastic element
	eta = 105e11*to.ones(mom_eq.n_elems)
	E1 = 10*ut.GPa*to.ones(mom_eq.n_elems)
	nu1 = 0.32*to.ones(mom_eq.n_elems)
	kelvin = sf.Viscoelastic(eta, E1, nu1, "kelvin")
```

### Block 15: Create creep

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
	# Create creep
	A = 1.9e-20*to.ones(mom_eq.n_elems)
	Q = 51600*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	creep_0 = sf.DislocationCreep(A, Q, n, "creep")
```

### Block 16: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
	# Create constitutive model
	mat.add_to_elastic(spring_0)
	mat.add_to_non_elastic(kelvin)
	mat.add_to_non_elastic(creep_0)
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

The temperature arrays define the reference and current thermal state used by temperature-dependent mechanical mechanisms.

```Python
	# Set initial temperature field
	T0_field = 298*to.ones(mom_eq.n_elems)
	mom_eq.set_T0(T0_field)
	mom_eq.set_T(T0_field)
```

### Block 20: Time settings for equilibrium stage

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
	# Time settings for equilibrium stage
	# tc_equilibrium = TimeControllerParabolic(final_time=10, initial_time=0.0, n_time_steps=50, time_unit="day")
	tc_equilibrium = sf.TimeController(dt=0.5, initial_time=0.0, final_time=10, time_unit="hour")
```

### Block 21: Boundary conditions

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Boundary conditions
	time_values = [0*ut.hour,  1*ut.hour]
	nt = len(time_values)
```

### Block 22: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_west = momBC.DirichletBC(boundary_name = "West",
							component = 0,
							values = [0.0, 0.0],
							time_values = [0.0, tc_equilibrium.t_final])
```

### Block 23: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_bottom = momBC.DirichletBC(boundary_name = "Bottom",
						  component = 2,
						  values = [0.0, 0.0],
						  time_values = [0.0, tc_equilibrium.t_final])
```

### Block 24: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_south = momBC.DirichletBC(boundary_name = "South",
						  component = 1,
						  values = [0.0, 0.0],
						  time_values = [0.0, tc_equilibrium.t_final])
```

### Block 25: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	side_burden = 10.0*ut.MPa
	bc_east = momBC.NeumannBC(boundary_name = "East",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0, tc_equilibrium.t_final],
						g = g_vec[2])
```

### Block 26: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_north = momBC.NeumannBC(boundary_name = "North",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0, tc_equilibrium.t_final],
						g = g_vec[2])
```

### Block 27: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	over_burden = 10.0*ut.MPa
	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [over_burden, over_burden],
						time_values = [0.0, tc_equilibrium.t_final],
						g = g_vec[2])
```

### Block 28: density = 0.082,

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	gas_density = 0.082
	p_gas = 10.0*ut.MPa
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						# density = 0.082,
						ref_pos = 430.0,
						values =      [p_gas, p_gas],
						time_values = [0.0,            tc_equilibrium.t_final],
						g = g_vec[2])
```

### Block 29: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_equilibrium = momBC.BcHandler(mom_eq)
	bc_equilibrium.add_boundary_condition(bc_west)
	bc_equilibrium.add_boundary_condition(bc_bottom)
	bc_equilibrium.add_boundary_condition(bc_south)
	bc_equilibrium.add_boundary_condition(bc_east)
	bc_equilibrium.add_boundary_condition(bc_north)
	bc_equilibrium.add_boundary_condition(bc_top)
	bc_equilibrium.add_boundary_condition(bc_cavern)
```

### Block 30: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_equilibrium)
```

### Block 31: Equilibrium output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Equilibrium output folder
	ouput_folder_equilibrium = os.path.join(output_folder, "equilibrium")
```

### Block 32: Print output folder

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(ouput_folder_equilibrium)
```

### Block 33: Create output handlers

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

### Block 34: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
	# Define simulator
	sim = sf.Simulator_M(mom_eq, tc_equilibrium, outputs, compute_elastic_response=True)
	sim.run()
```

### Block 35: Create Desai's viscoplastic model

The Desai viscoplastic element introduces yield-controlled inelastic deformation and hardening. It is added after the equilibrium stage when the example wants operational plasticity rather than only creep and viscoelasticity.

```Python




	# Create Desai's viscoplastic model
	mu_1 = 5.3665857009859815e-11*to.ones(mom_eq.n_elems)
	N_1 = 3.1*to.ones(mom_eq.n_elems)
	n = 3.0*to.ones(mom_eq.n_elems)
	a_1 = 1.965018496922832e-05*to.ones(mom_eq.n_elems)
	eta = 0.8275682807874163*to.ones(mom_eq.n_elems)
	beta_1 = 0.0048*to.ones(mom_eq.n_elems)
	beta = 0.995*to.ones(mom_eq.n_elems)
	m = -0.5*to.ones(mom_eq.n_elems)
	gamma = 0.095*to.ones(mom_eq.n_elems)
	alpha_0 = 0.0022*to.ones(mom_eq.n_elems)
	sigma_t = 5.0*to.ones(mom_eq.n_elems)
	desai = sf.ViscoplasticDesai(mu_1, N_1, a_1, eta, n, beta_1, beta, m, gamma, sigma_t, alpha_0, "desai")
```

### Block 36: Compute initial hardening parameter

The Desai model needs a consistent starting hardening variable. The script derives it from the stress state produced by the equilibrium stage so the operation stage does not begin with an artificial yield-function jump.

```Python
	# Compute initial hardening parameter
	stress_to = ut.numpy2torch(mom_eq.sig.x.array.reshape((mom_eq.n_elems, 3, 3)))
	desai.compute_initial_hardening(stress_to, Fvp_0=0.0)
```

### Block 37: Add viscoplastic element to constitutive model

These calls register history-dependent inelastic mechanisms with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
	# Add viscoplastic element to constitutive model
	mom_eq.mat.add_to_non_elastic(desai)
```

### Block 38: Time settings for operation stage

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python
	# Time settings for operation stage
	tc_operation = sf.TimeController(dt=0.1, initial_time=0.0, final_time=24, time_unit="hour")
```

### Block 39: Boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	# Boundary conditions
	bc_west = momBC.DirichletBC(boundary_name = "West",
							component = 0,
							values = [0.0, 0.0],
							time_values = [0.0, tc_operation.t_final])
```

### Block 40: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_bottom = momBC.DirichletBC(boundary_name = "Bottom",
						  component = 2,
						  values = [0.0, 0.0],
						  time_values = [0.0, tc_operation.t_final])
```

### Block 41: Continuation of the script

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
	bc_south = momBC.DirichletBC(boundary_name = "South",
						  component = 1,
						  values = [0.0, 0.0],
						  time_values = [0.0, tc_operation.t_final])
```

### Block 42: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_east = momBC.NeumannBC(boundary_name = "East",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0, tc_operation.t_final],
						g = g_vec[2])
```

### Block 43: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_north = momBC.NeumannBC(boundary_name = "North",
						direction = 2,
						density = salt_density,
						ref_pos = 660.0,
						values =      [side_burden, side_burden],
						time_values = [0.0, tc_operation.t_final],
						g = g_vec[2])
```

### Block 44: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_top = momBC.NeumannBC(boundary_name = "Top",
						direction = 2,
						density = 0.0,
						ref_pos = 0.0,
						values =      [over_burden, over_burden],
						time_values = [0.0, tc_operation.t_final],
						g = g_vec[2])
```

### Block 45: Continuation of the script

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
	bc_cavern = momBC.NeumannBC(boundary_name = "Cavern",
						direction = 2,
						density = gas_density,
						ref_pos = 430.0,
						values =      [10.0*ut.MPa, 7.0*ut.MPa, 7.0*ut.MPa, 10.0*ut.MPa, 10.0*ut.MPa],
						time_values = [0.0, 2.0*ut.hour, 14*ut.hour, 16*ut.hour, 24*ut.hour],
						g = g_vec[2])
```

### Block 46: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	bc_operation = momBC.BcHandler(mom_eq)
	bc_operation.add_boundary_condition(bc_west)
	bc_operation.add_boundary_condition(bc_bottom)
	bc_operation.add_boundary_condition(bc_south)
	bc_operation.add_boundary_condition(bc_east)
	bc_operation.add_boundary_condition(bc_north)
	bc_operation.add_boundary_condition(bc_top)
	bc_operation.add_boundary_condition(bc_cavern)
```

### Block 47: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
	# Set boundary conditions
	mom_eq.set_boundary_conditions(bc_operation)
```

### Block 48: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Define output folder
	output_folder_operation = os.path.join(output_folder, "operation")
```

### Block 49: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Print output folder
	if MPI.COMM_WORLD.rank == 0:
		print(output_folder_operation)
```

### Block 50: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
	# Create output handlers
	output_mom = sf.SaveFields(mom_eq)
	output_mom.set_output_folder(output_folder_operation)
	output_mom.add_output_field("u", "Displacement (m)")
	output_mom.add_output_field("eps_tot", "Total strain (-)")
	output_mom.add_output_field("eps_vp", "Viscoplastic strain (-)")
	output_mom.add_output_field("alpha", "Hardening parameter (-)")
	output_mom.add_output_field("Fvp", "Yield function (-)")
	output_mom.add_output_field("p_elems", "Mean stress (Pa)")
	output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
	outputs = [output_mom]
```

### Block 51: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
	# Define simulator
	sim = sf.Simulator_M(mom_eq, tc_operation, outputs, compute_elastic_response=False)
	sim.run()
```

### Block 52: run("P1")

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
	# run("P1")
	# run("P1P1")
	run("P1P1_Stab")
```

### Block 53: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python
if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: reads the cavern mesh and its boundary and region tags.
- `safeincave.LinearMomentum` and `safeincave.LinearMomentumMixed`: define the mechanical balance equation.
- `safeincave.Material`: stores density and the constitutive model used by the mechanics solver.
- `safeincave.Spring`, `safeincave.Viscoelastic`, `safeincave.DislocationCreep`, and `safeincave.ViscoplasticDesai`: define elastic, time-dependent, creep, and viscoplastic salt behavior.
- `safeincave.MomentumBC.DirichletBC` and `safeincave.MomentumBC.NeumannBC`: prescribe displacement constraints and pressure loads.
- `safeincave.TimeController`: advances the equilibrium and operation stages.
- `safeincave.SaveFields`: writes displacement, strain, stress, pressure, and von Mises stress fields.
- `safeincave.Simulator_M`: runs the mechanical simulation.

## Running the Example

From the repository root, run:

```bash
cd examples/mechanics/3_cavern
python main.py
```

The example uses `input_file.json` for configuration data and the mesh in `grids/cavern_regular`. Results are written to the local `output/` folder. To visualize the saved fields, run:

```bash
python plot_results.py
```

## Conclusion

This example demonstrates how to build a staged mechanical cavern simulation, initialize an equilibrium state, then continue with a different operation history and an enriched constitutive model. The same pattern can be used for pressure-cycling studies and for comparing different mechanical formulations.
