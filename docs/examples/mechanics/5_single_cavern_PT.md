# Example 5: Single Cavern With Prescribed Pressure and Temperature

## Introduction

This example couples the mechanical cavern model to a cavern boundary model with prescribed pressure and temperature histories. The focus is on mechanical deformation and stress driven by the specified cavern operation.

## Problem description

This section walks through `examples/mechanics/5_single_cavern_PT/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, cavern-fluid boundary models, PETSc linear solvers, field sampling utilities. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
from safeincave.Utils import day, GPa, create_field_elems
import safeincave.MomentumBC as momBC
import safeincave.CavernBC as caveBC
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


def main():
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cavern_regular")
    grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 4: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Define output folder
    output_folder = os.path.join("output", "case_0")
```

### Block 5: Define momentum equation

The mixed momentum formulation is selected when pressure-related stability or volumetric behavior matters. The `stab_scaling` parameter controls the pressure-stabilization contribution, while `theta` sets the time-integration convention used by inelastic updates.

```Python
    # Define momentum equation
    theta = 0.0
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

### Block 8: Set material density

The density vector is stored in the `Material` object. Mechanics uses it when building body forces, while thermal examples keep it together with heat capacity for the transient heat-storage term.

```Python
    # Set material density
    salt_density = 2200
    rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)
```

### Block 9: Constitutive model

The `Spring` element supplies the isotropic elastic stiffness matrix used as the backbone of the mechanical model. Other inelastic mechanisms evolve around this elastic response.

```Python
    # Constitutive model
    E0 = 102*GPa*to.ones(mom_eq.n_elems)
    nu0 = 0.3*to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")
```

### Block 10: Create creep

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
    # Create creep
    A = 1.9e-20*to.ones(mom_eq.n_elems)
    Q = 51600*to.ones(mom_eq.n_elems)
    n = 3.0*to.ones(mom_eq.n_elems)
    creep_0 = sf.DislocationCreep(A, Q, n, "creep")
```

### Block 11: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
    # Create constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_non_elastic(creep_0)
```

### Block 12: Set constitutive model

Attaching the material to the momentum equation initializes stiffness-related fields and makes the constitutive mechanisms available to the mechanical solver.

```Python
    # Set constitutive model
    mom_eq.set_material(mat)
```

### Block 13: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
    # Set body forces
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)
```

### Block 14: Set initial temperature field

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

### Block 15: Time settings for equilibrium stage

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
    # Time settings for equilibrium stage
    tc_eq = sf.TimeControllerParabolic(n_time_steps=100,
                                        initial_time=0.0,
                                        final_time=365,
                                        time_unit="day")
```

### Block 16: Boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
    # Boundary conditions
    bc_equilibrium = momBC.BcHandler(mom_eq)
```

### Block 17: Apply Dirichlet boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
    # Apply Dirichlet boundary conditions
    boundaries = [("West", 0),
                    ("East", 0),
                    ("South", 1),
                    ("North", 1),
                    ("Bottom", 2)]
    for b_name, component in boundaries:
        bc = momBC.DirichletBC(boundary_name=b_name, component=component, values=[0.0, 0.0], time_values=[0.0, tc_eq.t_final])
        bc_equilibrium.add_boundary_condition(bc)
```

### Block 18: Apply overburden

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
    # Apply overburden
    overburden = 10*sf.Utils.MPa
    bc_top = momBC.NeumannBC(boundary_name = "Top",
                        direction = 2,
                        density = 0.0,
                        ref_pos = 0.0,
                        values = [overburden, overburden],
                        time_values = [0*day,  tc_eq.t_final],
                        g = g_vec[2])
    bc_equilibrium.add_boundary_condition(bc_top)
```

### Block 19: Calculate lithostatic pressure at the cavern's roof

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

    # Calculate lithostatic pressure at the cavern's roof
    hanging_wall = 430
    p_roof = overburden + salt_density*abs(g)*hanging_wall
```

### Block 20: Define cavern conditions

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python

    # Define cavern conditions
    cavern_handler = caveBC.CavernHandler()
    cave_1 = caveBC.Cavern_PT(
                            grid = grid,
                            cavern_name = "Cavern",
                            fluid = "Hydrogen",
                            sym_scale = 4,
                            reference_point = [0.0, 0.0, hanging_wall],
                            P_values = [0.8*p_roof, 0.8*p_roof],
                            T_values = [303, 303],
                            time_values = [0*day,  tc_eq.t_final],
                            ref_pos = hanging_wall,
                            direction = 2,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_1)
    cavern_handler.set_output_folder(output_folder)
```

### Block 21: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_equilibrium)
```

### Block 22: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Print output folder
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder)
        sys.stdout.flush()
```

### Block 23: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs = [output_mom]
```

### Block 24: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
    # Define simulator
    sim = sf.Simulator_M(mom_eq, tc_eq, outputs, cavern_handler, True)
    sim.run()
```

### Block 25: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python




if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.LinearMomentumMixed`: solves the stabilized mixed mechanical formulation.
- `safeincave.CavernBC.Cavern_PT`: applies prescribed cavern pressure and temperature histories.
- `safeincave.CavernBC.CavernHandler`: connects the cavern model to the simulator.
- `safeincave.Spring` and `safeincave.DislocationCreep`: define elastic and creep behavior.
- `safeincave.TimeControllerParabolic`: defines the operation time grid.
- `safeincave.MomentumBC`: sets displacement constraints and overburden loading.
- `safeincave.Simulator_M`: runs the mechanical model with cavern coupling.

## Running the Example

From the repository root, run:

```bash
cd examples/mechanics/5_single_cavern_PT
python main.py
```

Results are written to `output/`. Plot cavern histories with:

```bash
python plot_cavern_data.py
```

## Conclusion

This example shows how to drive a mechanical cavern simulation with a prescribed pressure-temperature operation. It is useful for mechanical assessment when the cavern fluid history is known from measurements or another simulator.
