# Example 7: Multi-Cavern Mechanics

## Introduction

This example models a multicavern system with different cavern operation types. One cavern uses methane mass flux, one uses prescribed hydrogen pressure and temperature, and one uses water mass flux.

## Problem description

This section walks through `examples/mechanics/7_multicavern/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

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
    grid_path = os.path.join("..", "..", "..", "grids", "multicavern")
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

### Block 15: Define time controler

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python
    # Define time controler
    t_0 = 0.0
    t_final = 300
    dt = t_final/600
    tc = sf.TimeController(dt=dt, initial_time=t_0, final_time=t_final, time_unit="day")
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
    boundaries = [("West", 0), ("East", 0), ("South", 1), ("North", 1), ("Bottom", 2)]
    for b_name, component in boundaries:
        bc = momBC.DirichletBC(boundary_name=b_name, component=component, values=[0.0, 0.0], time_values=[0.0, tc.t_final])
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
                        time_values = [0*day,  tc.t_final],
                        g = g_vec[2])
    bc_equilibrium.add_boundary_condition(bc_top)
```

### Block 19: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_equilibrium)
```

### Block 20: Calculate lithostatic pressure at cavern's midpoint

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

    # Calculate lithostatic pressure at cavern's midpoint
    z_roof = 420
    z_floor = 200
    z_ground = 660
    z_mid = (z_roof + z_floor) / 2
    p_mid = overburden + salt_density*abs(g)*(z_ground - z_mid)
    p_roof = overburden + salt_density*abs(g)*z_roof
    print(0.2*p_roof/sf.Utils.MPa, 0.8*p_roof/sf.Utils.MPa)
```

### Block 21: Read flow rate values

The operation schedule is read from JSON so pressure, temperature, or mass-flow histories can be changed without editing solver setup code. The script then passes those arrays directly into the cavern model.

```Python

    # Read flow rate values
    data_caverns = sf.Utils.read_json("input_cavern_data.json")
```

### Block 22: Define caverns

`CavernHandler` collects all cavern boundary models and gives the simulator one interface for volume updates, heat-transfer integration, pressure loads, and output histories. This is especially important for multicavern examples.

```Python
    # Define caverns
    cavern_handler = caveBC.CavernHandler()
```

### Block 23: Continuation of the script

`Cavern_MassFlux` is selected when injection or withdrawal rates drive the cavern state. It combines the mass-flow schedule, heat exchange, current volume, and CoolProp thermodynamics to update pressure, temperature, density, and mass.

```Python
    cave_methane = caveBC.Cavern_MassFlux(
                            grid = grid,
                            cavern_name = "Cavern_full",
                            fluid = data_caverns["Cavern_full"]["fluid"],
                            sym_scale = 1,
                            reference_point = [400, 400, 330],
                            P_init = data_caverns["Cavern_full"]["P_init"],
                            T_init = data_caverns["Cavern_full"]["T_init"],
                            T_in = 300.0,
                            Mflux_values = data_caverns["Cavern_full"]["flow"],
                            time_values = data_caverns["Cavern_full"]["time"],
                            direction = 2,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_methane)
```

### Block 24: Continuation of the script

`Cavern_PT` is used when pressure and temperature histories are known externally. It still computes cavern volume and records data, but it does not solve a mass-balance problem for pressure.

```Python
    cave_hydrogen = caveBC.Cavern_PT(
                            grid = grid,
                            cavern_name = "Cavern_half",
                            fluid = data_caverns["Cavern_half"]["fluid"],
                            sym_scale = 2,
                            reference_point = [800, 400, 330],
                            P_values = data_caverns["Cavern_half"]["P_hist"],
                            T_values = data_caverns["Cavern_half"]["T_hist"],
                            time_values = data_caverns["Cavern_half"]["time"],
                            ref_pos = z_roof,
                            direction = 2,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_hydrogen)
```

### Block 25: Continuation of the script

`Cavern_MassFlux` is selected when injection or withdrawal rates drive the cavern state. It combines the mass-flow schedule, heat exchange, current volume, and CoolProp thermodynamics to update pressure, temperature, density, and mass.

```Python
    cave_water = caveBC.Cavern_MassFlux(
                            grid = grid,
                            cavern_name = "Cavern_quarter",
                            fluid = data_caverns["Cavern_quarter"]["fluid"],
                            sym_scale = 4,
                            reference_point = [0.0, 0.0, 330],
                            P_init = data_caverns["Cavern_quarter"]["P_init"],
                            T_init = data_caverns["Cavern_quarter"]["T_init"],
                            T_in = 0.0,
                            Mflux_values = data_caverns["Cavern_quarter"]["flow"],
                            time_values = data_caverns["Cavern_quarter"]["time"],
                            direction = 2,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_water)
```

### Block 26: Continuation of the script

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    cavern_handler.set_output_folder(output_folder)
```

### Block 27: Print output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Print output folder
    if MPI.COMM_WORLD.rank == 0:
        print(output_folder)
        sys.stdout.flush()
```

### Block 28: Create output handlers

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

### Block 29: Define simulator

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
    # Define simulator
    sim = sf.Simulator_M(mom_eq, tc, outputs, cavern_handler, True)
    sim.run()
```

### Block 30: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python




if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: loads the multicavern mesh.
- `safeincave.Utils.read_json`: reads `input_cavern_data.json`.
- `safeincave.CavernBC.Cavern_MassFlux`: represents caverns controlled by mass exchange.
- `safeincave.CavernBC.Cavern_PT`: represents a cavern controlled by pressure and temperature histories.
- `safeincave.CavernBC.CavernHandler`: groups all cavern models.
- `safeincave.LinearMomentumMixed`: solves the stabilized mixed mechanics problem.
- `safeincave.Spring` and `safeincave.DislocationCreep`: define the salt mechanical behavior.
- `safeincave.Simulator_M`: runs the multicavern mechanical simulation.

## Running the Example

From the repository root, run:

```bash
cd examples/mechanics/7_multicavern
python main.py
```

The helper script can regenerate operation data:

```bash
python create_cavern_operations.py
```

After the run, plot cavern histories with:

```bash
python plot_cavern_data.py
```

## Conclusion

This example shows how multiple cavern boundary models can be assembled into one mechanical simulation. It is a starting point for comparing interactions between nearby storage caverns.
