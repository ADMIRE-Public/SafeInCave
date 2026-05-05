# Example 5: Methane Cavern With Mass Flux

## Introduction

This example simulates a methane storage cavern with prescribed mass-flux data. The cavern thermodynamic state is updated from injection and withdrawal data, and the resulting pressure and temperature act on the coupled rock model.

## Problem description

This section walks through `examples/thermomechanics/5_cavern_methane/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

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
    output_folder = os.path.join("output", "case_1")
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

### Block 7: Define heat diffusion equation

`HeatDiffusion` allocates the scalar temperature space, cell-wise material fields, UFL measures, and temperature storage used by the thermal solve. The example creates it before assigning material properties because `set_material` copies those properties into DG0 finite-element fields.

```Python
    # Define heat diffusion equation
    heat_eq = sf.HeatDiffusion(grid)
```

### Block 8: Define solver

The thermal linear system is solved with PETSc. Conjugate gradients match the symmetric diffusion operator, and the additive Schwarz preconditioner keeps the setup usable in both serial and parallel runs.

```Python
    # Define solver
    solver_heat = PETSc.KSP().create(grid.mesh.comm)
    solver_heat.setType("cg")
    solver_heat.getPC().setType("asm")
    solver_heat.setTolerances(rtol=1e-12, max_it=100)
    heat_eq.set_solver(solver_heat)
```

### Block 9: Define material properties

`Material` is the central container for element-wise properties and constitutive mechanisms. The equation object later asks this container for stiffness, density, creep operators, thermal expansion, and thermal transport data.

```Python
    # Define material properties
    mat = sf.Material(mom_eq.n_elems)
```

### Block 10: Set material density

The density vector is stored in the `Material` object. Mechanics uses it when building body forces, while thermal examples keep it together with heat capacity for the transient heat-storage term.

```Python
    # Set material density
    salt_density = 2200
    rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)
```

### Block 11: Constitutive model

The `Spring` element supplies the isotropic elastic stiffness matrix used as the backbone of the mechanical model. Other inelastic mechanisms evolve around this elastic response.

```Python
    # Constitutive model
    E0 = 102*GPa*to.ones(mom_eq.n_elems, dtype=to.float64)
    nu0 = 0.3*to.ones(mom_eq.n_elems, dtype=to.float64)
    spring_0 = sf.Spring(E0, nu0, "spring")
```

### Block 12: Create creep

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
    # Create creep
    A = 1.9e-20*to.ones(mom_eq.n_elems, dtype=to.float64)
    Q = 51600*to.ones(mom_eq.n_elems, dtype=to.float64)
    n = 3.0*to.ones(mom_eq.n_elems, dtype=to.float64)
    creep_ds = sf.DislocationCreep(A, Q, n, "ds_creep")
```

### Block 13: Create pressure solution creep

`PressureSolutionCreep` represents another salt-creep pathway, controlled by stress, temperature, activation energy, and grain-size parameter `d`. It is combined with dislocation creep through the `Material` non-elastic element list.

```Python
    # Create pressure solution creep
    A = 1.29e-29*to.ones(mom_eq.n_elems, dtype=to.float64)
    Q = 13184*to.ones(mom_eq.n_elems, dtype=to.float64)
    d = 0.01*to.ones(mom_eq.n_elems, dtype=to.float64)
    creep_ps = sf.PressureSolutionCreep(A, d, Q, "ps_creep")
```

### Block 14: Thermo-elastic element

The thermoelastic element converts temperature changes into thermal strain, `alpha * deltaT * I`. Adding it to the material is what lets heat diffusion produce mechanical stress and deformation in coupled examples.

```Python
    # Thermo-elastic element
    # alpha = 44e-7*to.ones(mom_eq.n_elems)
    alpha = 120e-6*to.ones(mom_eq.n_elems)
    thermo = sf.Thermoelastic(alpha, "thermo")
```

### Block 15: Set specific heat capacity

Specific heat capacity controls thermal inertia. SafeInCave stores it on the material and copies it into a DG0 field inside `HeatDiffusion`, where it multiplies the transient temperature term.

```Python
    # Set specific heat capacity
    cp = 850*to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_specific_heat_capacity(cp)
```

### Block 16: Set thermal conductivity

Thermal conductivity sets the strength of the diffusion operator. Defining it per element keeps the example ready for heterogeneous thermal properties even when the current values are uniform.

```Python
    # Set thermal conductivity
    k = 7*to.ones(heat_eq.n_elems, dtype=to.float64)
    mat.set_thermal_conductivity(k)
```

### Block 17: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms, thermal strain coupling with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
    # Create constitutive model
    mat.add_to_elastic(spring_0)
    mat.add_to_thermoelastic(thermo)
    mat.add_to_non_elastic(creep_ds)
    mat.add_to_non_elastic(creep_ps)
```

### Block 18: Set material properties to governing equations

The same material object is attached to both equations so thermal transport, thermal expansion, and mechanical constitutive response use consistent element-wise data.

```Python
    # Set material properties to governing equations
    mom_eq.set_material(mat)
    heat_eq.set_material(mat)
```

### Block 19: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
    # Set body forces
    g = -9.81
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)
```

### Block 20: Set initial temperature field

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

### Block 21: Time settings

This fixed-step controller sets the operational cadence explicitly. It is useful when the imposed pressure, temperature, or mass-flow schedule should be sampled with a predictable step size.

```Python
    # Time settings
    time_ctrl = sf.TimeController(dt=1.0, initial_time=0.0, final_time=365, time_unit="day")
    time_values = [time_ctrl.t_initial, time_ctrl.t_final]
    nt = len(time_values)
```

### Block 22: Boundary conditions for momentum balance equation

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
    # Boundary conditions for momentum balance equation
    bc_mom = momBC.BcHandler(mom_eq)
```

### Block 23: Apply Dirichlet boundary conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
    # Apply Dirichlet boundary conditions
    boundaries = [("West", 0),
                    ("East", 0),
                    ("South", 1),
                    ("North", 1),
                    ("Bottom", 2)]
    for b_name, component in boundaries:
        bc = momBC.DirichletBC(boundary_name=b_name, component=component, values=nt*[0.0], time_values=time_values)
        bc_mom.add_boundary_condition(bc)
```

### Block 24: Apply overburden

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

### Block 25: Set boundary conditions for momentum balance equation

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
    # Set boundary conditions for momentum balance equation
    mom_eq.set_boundary_conditions(bc_mom)
```

### Block 26: Set initial temperature

The initial temperature initializes both current and previous thermal states inside `HeatDiffusion`. Without this step the transient solve would not have a physically meaningful starting field.

```Python
    # Set initial temperature
    T0_field_nodes = create_field_nodes(grid, T_field_fun)
    heat_eq.set_initial_T(T0_field_nodes)
```

### Block 27: Boundary conditions for heat diffusion equation

This Dirichlet condition fixes temperature on a named boundary. The values are time-interpolated by the heat `BcHandler`, which turns the user schedule into DOLFINx essential boundary conditions at each step.

```Python
    # Boundary conditions for heat diffusion equation
    bc_heat = heatBC.BcHandler(heat_eq)
    bc_top = heatBC.DirichletBC("Top", nt*[T_top], time_values)
    bc_bottom = heatBC.NeumannBC("Bottom", nt*[dTdZ], time_values)
    bc_west = heatBC.NeumannBC("West", nt*[0.0], time_values)
    bc_east = heatBC.NeumannBC("East", nt*[0.0], time_values)
    bc_north = heatBC.NeumannBC("North", nt*[0.0], time_values)
    bc_south = heatBC.NeumannBC("South", nt*[0.0], time_values)
    bc_heat.add_boundary_condition(bc_top)
    bc_heat.add_boundary_condition(bc_bottom)
    bc_heat.add_boundary_condition(bc_west)
    bc_heat.add_boundary_condition(bc_east)
    bc_heat.add_boundary_condition(bc_north)
    bc_heat.add_boundary_condition(bc_south)
    heat_eq.set_boundary_conditions(bc_heat)
```

### Block 28: Calculate lithostatic pressure at the cavern's roof

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

    # Calculate lithostatic pressure at the cavern's roof
    hanging_wall = 430
    p_roof = overburden + salt_density*abs(g)*hanging_wall
    print(0.5*p_roof)
```

### Block 29: Read flow rate values

The operation schedule is read from JSON so pressure, temperature, or mass-flow histories can be changed without editing solver setup code. The script then passes those arrays directly into the cavern model.

```Python
    # Read flow rate values
    data_cavern = sf.Utils.read_json("input_massflux.json")
```

### Block 30: Define cavern conditions

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Define cavern conditions
    cavern_handler = caveBC.CavernHandler()
    cave_0 = caveBC.Cavern_MassFlux(
                            grid = grid,
                            cavern_name = "Cavern",
                            fluid = "Methane",
                            sym_scale = 4,
                            reference_point = [0.0, 0.0, hanging_wall],
                            P_init = 0.5*p_roof,
                            T_init = T_top,
                            T_in = T_top,
                            Mflux_values = data_cavern["flow"],
                            time_values = data_cavern["time"],
                            direction = 2,
                            h_conv = 5.0,
                            g = g_vec[2])
    cavern_handler.add_cavern(cave_0)
    cavern_handler.set_output_folder(output_folder)
```

### Block 31: Create output handlers

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

### Block 32: Continuation of the script

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    output_heat = sf.SaveFields(heat_eq)
    output_heat.set_output_folder(output_folder)
    output_heat.add_output_field("T", "Temperature (K)")
```

### Block 33: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
    outputs = [output_mom, output_heat]
```

### Block 34: Define simulator

`Simulator_Full` runs the most complete loop: heat diffusion, cavern heat exchange, cavern thermodynamics, cavern pressure loads, mechanical solves, and inelastic internal-variable updates are advanced together.

```Python
    # Define simulator
    sim = sf.Simulator_Full(mom_eq, heat_eq, time_ctrl, outputs, cavern_handler, True)
    # sim = sf.Simulator_M(mom_eq, time_ctrl, outputs, cavern_handler, True)
    sim.run()
```

### Block 35: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python




if __name__ == '__main__':
	main()
```

## Main Functions

- `safeincave.Utils.read_json`: reads `input_massflux.json`.
- `safeincave.CavernBC.Cavern_MassFlux`: updates methane pressure and temperature from mass flux, initial state, and inlet temperature.
- `safeincave.CavernBC.CavernHandler`: collects the cavern boundary model.
- `safeincave.LinearMomentumMixed` and `safeincave.HeatDiffusion`: solve the coupled rock mechanics and heat equations.
- `safeincave.Spring`, `safeincave.DislocationCreep`, `safeincave.PressureSolutionCreep`, and `safeincave.Thermoelastic`: define the salt response.
- `safeincave.Simulator_Full`: runs the coupled cavern-fluid, thermal, and mechanical simulation.

## Running the Example

From the repository root, run:

```bash
cd examples/thermomechanics/5_cavern_methane
python main.py
```

If the input schedule needs to be regenerated, use:

```bash
python build_massflow_data.py
```

After the simulation, plot the cavern data with:

```bash
python plot_cavern_data.py
```

## Conclusion

This example demonstrates operational methane storage with mass-flux-controlled cavern thermodynamics. It is a practical template for injection-withdrawal schedules coupled to rock deformation and heat transfer.
