# Example 7: Multi-Cavern Thermomechanics

## Introduction

This example simulates a system with multiple caverns using different fluids and operation models. Methane and water caverns use mass-flux histories, while the hydrogen cavern uses prescribed pressure and temperature data.

## Problem description

This section walks through `examples/thermomechanics/7_multicavern/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, thermal boundary conditions, cavern-fluid boundary models, PETSc linear solvers, field sampling utilities. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
from safeincave.Utils import day, GPa, create_field_elems, create_field_nodes
import safeincave.MomentumBC as momBC
import safeincave.HeatBC as heatBC
import safeincave.CavernBC as caveBC
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
    grid_path = os.path.join("..", "..", "..", "grids", "multicavern")
    grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 3: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Define output folder
    output_folder = os.path.join("output", "case_2")
```

### Block 4: Define momentum equation

The mixed momentum formulation is selected when pressure-related stability or volumetric behavior matters. The `stab_scaling` parameter controls the pressure-stabilization contribution, while `theta` sets the time-integration convention used by inelastic updates.

```Python

    # Define momentum equation
    mom_eq = sf.LinearMomentumMixed(grid, theta=1.0, stab_scaling=1.0)
    mom_solver = PETSc.KSP().create(grid.mesh.comm)
    mom_solver.setType("gmres")
    mom_solver.getPC().setType("asm")
    mom_solver.setTolerances(rtol=1e-12, max_it=100)
    mom_eq.set_solver(mom_solver)
```

### Block 5: Define heat diffusion equation

`HeatDiffusion` allocates the scalar temperature space, cell-wise material fields, UFL measures, and temperature storage used by the thermal solve. The example creates it before assigning material properties because `set_material` copies those properties into DG0 finite-element fields.

```Python

    # Define heat diffusion equation
    heat_eq = sf.HeatDiffusion(grid)
    solver_heat = PETSc.KSP().create(grid.mesh.comm)
    solver_heat.setType("cg")
    solver_heat.getPC().setType("asm")
    solver_heat.setTolerances(rtol=1e-12, max_it=100)
    heat_eq.set_solver(solver_heat)
```

### Block 6: Define material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python


    # Define material properties
    mat = sf.Material(mom_eq.n_elems)
```

### Block 7: Set material density

The density vector is stored in the `Material` object. Mechanics uses it when building body forces, while thermal examples keep it together with heat capacity for the transient heat-storage term.

```Python
    # Set material density
    salt_density = 2200
    rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)
```

### Block 8: Constitutive model

The `Spring` element supplies the isotropic elastic stiffness matrix used as the backbone of the mechanical model. Other inelastic mechanisms evolve around this elastic response.

```Python
    # Constitutive model
    E0 = 102*GPa*to.ones(mom_eq.n_elems)
    nu0 = 0.3*to.ones(mom_eq.n_elems)
    spring_0 = sf.Spring(E0, nu0, "spring")
    mat.add_to_elastic(spring_0)
```

### Block 9: Create dislocation creep

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
    # Create dislocation creep
    A = 1.9e-20*to.ones(mom_eq.n_elems)
    Q = 51600*to.ones(mom_eq.n_elems)
    n = 3.0*to.ones(mom_eq.n_elems)
    creep_ds = sf.DislocationCreep(A, Q, n, "creep")
    mat.add_to_non_elastic(creep_ds)
```

### Block 10: Create pressure solution creep

`PressureSolutionCreep` represents another salt-creep pathway, controlled by stress, temperature, activation energy, and grain-size parameter `d`. It is combined with dislocation creep through the `Material` non-elastic element list.

```Python
    # Create pressure solution creep
    A = 1.29e-29*to.ones(mom_eq.n_elems, dtype=to.float64)
    Q = 13184*to.ones(mom_eq.n_elems, dtype=to.float64)
    d = 0.01*to.ones(mom_eq.n_elems, dtype=to.float64)
    creep_ps = sf.PressureSolutionCreep(A, d, Q, "ps_creep")
    mat.add_to_non_elastic(creep_ps)
```

### Block 11: Create Thermo-elastic element

The thermoelastic element converts temperature changes into thermal strain, `alpha * deltaT * I`. Adding it to the material is what lets heat diffusion produce mechanical stress and deformation in coupled examples.

```Python
    # Create Thermo-elastic element
    alpha = 120e-6*to.ones(mom_eq.n_elems)
    thermo = sf.Thermoelastic(alpha, "thermo")
    mat.add_to_thermoelastic(thermo)
```

### Block 12: Set specific heat capacity

Specific heat capacity controls thermal inertia. SafeInCave stores it on the material and copies it into a DG0 field inside `HeatDiffusion`, where it multiplies the transient temperature term.

```Python

    # Set specific heat capacity
    cp = 850*to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_specific_heat_capacity(cp)
```

### Block 13: Set thermal conductivity

Thermal conductivity sets the strength of the diffusion operator. Defining it per element keeps the example ready for heterogeneous thermal properties even when the current values are uniform.

```Python
    # Set thermal conductivity
    k = 7*to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_thermal_conductivity(k)
```

### Block 14: Set material properties to governing equations

The same material object is attached to both equations so thermal transport, thermal expansion, and mechanical constitutive response use consistent element-wise data.

```Python
    # Set material properties to governing equations
    mom_eq.set_material(mat)
    heat_eq.set_material(mat)
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

The geothermal temperature is sampled twice because the two solvers store temperature differently: mechanics consumes element-centered values for constitutive updates, while heat diffusion evolves nodal P1 temperatures.

```Python
    # Set initial temperature field
    km = 1000
    dTdZ = 27/km
    T_top = 273 + 20
    T_field_fun = lambda x,y,z: T_top + dTdZ*(660 - z)
    T0_field_elems = create_field_elems(grid, T_field_fun)
    T0_field_nodes = create_field_nodes(grid, T_field_fun)
    mom_eq.set_T0(T0_field_elems)
    mom_eq.set_T(T0_field_elems)
    heat_eq.set_initial_T(T0_field_nodes)
```

### Block 17: Define time controler

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python
    # Define time controler
    time_ctrl = sf.TimeController(dt=0.1, initial_time=1.0, final_time=365, time_unit="day")
    time_values = [time_ctrl.t_initial, time_ctrl.t_final]
    nt = len(time_values)
```

### Block 18: Boundary conditions for momemtum equation

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
    # Boundary conditions for momemtum equation
    bc_mom = momBC.BcHandler(mom_eq)
```

### Block 19: Apply Dirichlet boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
    # Apply Dirichlet boundary conditions
    boundaries = [("West", 0), ("East", 0), ("South", 1), ("North", 1), ("Bottom", 2)]
    for b_name, component in boundaries:
        bc = momBC.DirichletBC(boundary_name=b_name, component=component, values=nt*[0.0], time_values=time_values)
        bc_mom.add_boundary_condition(bc)
```

### Block 20: Apply overburden

This Neumann condition applies a pressure or traction history. The base value is interpolated in time, and the optional density/reference-position terms add hydrostatic variation along the selected coordinate direction.

```Python
    # Apply overburden
    overburden = 10*sf.Utils.MPa
    bc_top = momBC.NeumannBC(boundary_name = "Top",
                        direction = 2,
                        density = 0.0,
                        ref_pos = 0.0,
                        values = nt*[overburden],
                        time_values = time_values,
                        g = g_vec[2])
    bc_mom.add_boundary_condition(bc_top)
```

### Block 21: Set boundary conditions

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python

    # Set boundary conditions
    mom_eq.set_boundary_conditions(bc_mom)
```

### Block 22: Calculate lithostatic pressure at cavern's midpoint

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

### Block 23: Boundary conditions for heat diffusion equation

This Dirichlet condition fixes temperature on a named boundary. The values are time-interpolated by the heat `BcHandler`, which turns the user schedule into DOLFINx essential boundary conditions at each step.

```Python


    # Boundary conditions for heat diffusion equation
    bc_heat = heatBC.BcHandler(heat_eq)
    bc_top = heatBC.DirichletBC("Top",      [T_top, T_top], time_values)
    bc_bottom = heatBC.NeumannBC("Bottom",  [dTdZ, dTdZ],   time_values)
    bc_west = heatBC.NeumannBC("West",      [0.0, 0.0],     time_values)
    bc_east = heatBC.NeumannBC("East",      [0.0, 0.0],     time_values)
    bc_north = heatBC.NeumannBC("North",    [0.0, 0.0],     time_values)
    bc_south = heatBC.NeumannBC("South",    [0.0, 0.0],     time_values)
    bc_heat.add_boundary_condition(bc_top)
    bc_heat.add_boundary_condition(bc_bottom)
    bc_heat.add_boundary_condition(bc_west)
    bc_heat.add_boundary_condition(bc_east)
    bc_heat.add_boundary_condition(bc_north)
    bc_heat.add_boundary_condition(bc_south)
    heat_eq.set_boundary_conditions(bc_heat)
```

### Block 24: Calculate lithostatic pressure at cavern level

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

    # Calculate lithostatic pressure at cavern level
    z_roof = 420
    z_floor = 200
    z_ground = 660
    z_mid = (z_roof + z_floor) / 2
    p_mid = overburden + salt_density*abs(g)*(z_ground - z_mid)
    p_roof = overburden + salt_density*abs(g)*z_roof
    print(0.2*p_roof/sf.Utils.MPa, 0.8*p_roof/sf.Utils.MPa)
```

### Block 25: Define caverns

`CavernHandler` collects all cavern boundary models and gives the simulator one interface for volume updates, heat-transfer integration, pressure loads, and output histories. This is especially important for multicavern examples.

```Python


    # Define caverns
    cavern_handler = caveBC.CavernHandler()
```

### Block 26: Continuation of the script

`Cavern_MassFlux` is selected when injection or withdrawal rates drive the cavern state. It combines the mass-flow schedule, heat exchange, current volume, and CoolProp thermodynamics to update pressure, temperature, density, and mass.

```Python
    data_methane = sf.Utils.read_json("input_methane.json")
    cave_methane = caveBC.Cavern_MassFlux(
                            grid = grid,
                            cavern_name = "Cavern_full",
                            fluid = "Methane",
                            sym_scale = 1,
                            reference_point = [0.0, 0.0, z_roof],
                            P_init = 0.5*p_roof,
                            T_init = T_top,
                            T_in = T_top,
                            Mflux_values = data_methane["flow"],
                            time_values = data_methane["time"],
                            direction = 2,
                            h_conv = 5.0,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_methane)
```

### Block 27: Continuation of the script

`Cavern_PT` is used when pressure and temperature histories are known externally. It still computes cavern volume and records data, but it does not solve a mass-balance problem for pressure.

```Python
    data_hydrogen = sf.Utils.read_json("input_hydrogen.json")
    cave_hydrogen = caveBC.Cavern_PT(
                            grid = grid,
                            cavern_name = "Cavern_half",
                            fluid = "Hydrogen",
                            sym_scale = 2,
                            reference_point = [0.0, 0.0, z_roof],
                            P_values = data_hydrogen["pressure"],
                            T_values = data_hydrogen["temperature"],
                            time_values = data_hydrogen["time"],
                            ref_pos = z_roof,
                            direction = 2,
                            h_conv = 5.0,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_hydrogen)
```

### Block 28: Continuation of the script

`Cavern_MassFlux` is selected when injection or withdrawal rates drive the cavern state. It combines the mass-flow schedule, heat exchange, current volume, and CoolProp thermodynamics to update pressure, temperature, density, and mass.

```Python
    cave_water = caveBC.Cavern_MassFlux(
                            grid = grid,
                            cavern_name = "Cavern_quarter",
                            fluid = "Water",
                            sym_scale = 4,
                            reference_point = [0.0, 0.0, z_roof],
                            P_init = 0.5*p_roof,
                            T_init = T_top,
                            T_in = 0.0,
                            Mflux_values = [0.0, 0.0],
                            time_values = time_values,
                            direction = 2,
                            h_conv = 5.0,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_water)
```

### Block 29: Continuation of the script

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    cavern_handler.set_output_folder(output_folder)
```

### Block 30: Create output handlers

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python

    # Create output handlers
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(output_folder)
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
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

### Block 33: Define simulator

`Simulator_Full` runs the most complete loop: heat diffusion, cavern heat exchange, cavern thermodynamics, cavern pressure loads, mechanical solves, and inelastic internal-variable updates are advanced together.

```Python
    # Define simulator
    sim = sf.Simulator_Full(mom_eq, heat_eq, time_ctrl, outputs, cavern_handler, True)
    sim.run()
```

### Block 34: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python




if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: loads the multicavern mesh.
- `safeincave.Utils.read_json`: reads cavern geometry and operation input files.
- `safeincave.CavernBC.Cavern_MassFlux`: models caverns controlled by injection and withdrawal rates.
- `safeincave.CavernBC.Cavern_PT`: models a cavern controlled by pressure and temperature histories.
- `safeincave.CavernBC.CavernHandler`: manages all cavern boundary models.
- `safeincave.LinearMomentumMixed` and `safeincave.HeatDiffusion`: solve the coupled rock problem.
- `safeincave.Simulator_Full`: advances all cavern, thermal, and mechanical states together.

## Running the Example

From the repository root, run:

```bash
cd examples/thermomechanics/7_multicavern
python main.py
```

The example reads `input_cavern_data.json`, `input_methane.json`, and `input_hydrogen.json`. Results are saved under `output/`. Plot cavern histories with:

```bash
python plot_cavern_data.py
```

## Conclusion

This example demonstrates how SafeInCave can couple several cavern boundary models in one simulation. The same workflow can be extended to field-scale storage studies with different fluids, schedules, and cavern geometries.
