# Example 8: Surface Settlement

## Introduction

This example computes deformation and stress for a geological model intended for settlement analysis. It uses elastic and dislocation-creep behavior and saves displacement fields that can be used to evaluate surface movement.

## Problem description

This section walks through `examples/mechanics/8_settlement/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The imports reveal the physics activated by this example: mechanical boundary conditions, PETSc linear solvers, field sampling utilities. Loading these modules up front keeps the script explicit about which SafeInCave subsystems will be wired together.

```Python
import safeincave as sf
from safeincave.Utils import hour, day, GPa, MPa, create_field_elems
import safeincave.MomentumBC as momBC
from mpi4py import MPI
from petsc4py import PETSc
import torch as to
import os
import sys
```

### Block 2: Function `calculate_lithostatic_pressure`

This helper expresses the geostatic load as the weight of overlying layers. In a settlement-style calculation, that pressure is the natural scale for boundary loads or cavern-pressure ratios.

```Python
def calculate_lithostatic_pressure(salt_density, cap_density, ovb_density):
    '''
    This function calculates the lithostatic pressure at cavern's roof
```

### Block 3: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
      z ^                            z ^
        |                              |
        | __________  _ _ z0=0 _ _ _ _ |
        |           |                  |\
        |           |                  | \
        |  OVB      |                  |  \
        |           |                  |   \
        |           |                  |    \
        | __________| _ _z1=-670.56 m _|_ _ _\
        |           |                  |      `
        |  CAP      |                  |        `
        | __________| _ _ z2=-914.4 m _| _ _ _ _ _`
        |           |                  |            `
        |  SALT     |                  |             \
        |           |                  |              \
        | _ _ _ _ _ | _ z3=-1158.24 m _|_ _ _ _ _ _ _ _\
        |  \        |                  |               :\
        |   \       |                  |               : \
        |    |      |                  |               :  \
        |    |      |                  |               :   \
        |    |      |                  +---------------:------------> Lithostatic pressure
        |    |      |                                  P
        |   /       |
        | _/        |
        |           |
        |           |
        |___________|
    '''
    z0 = 0
    z1 = -670.56
    z2 = -914.4
    z3 = -1158.24
```

### Block 4: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
    H_ovb = z0 - z1
    H_cap = z1 - z2
    H_sal = z2 - z3
```

### Block 5: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
    g = 9.81
```

### Block 6: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
    P = (ovb_density*H_ovb + cap_density*H_cap + salt_density*H_sal) * g
    return P
```

### Block 7: Read grid

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python


def run(case_name):
    # Read grid
    grid_path = os.path.join("..", "..", "..", "grids", "cube")
    grid = sf.GridHandlerGMSH("geom", grid_path)
```

### Block 8: Define output folder

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Define output folder
    output_folder = os.path.join("output", f"case_{case_name}")
```

### Block 9: Define momentum equation

The displacement-based momentum equation is enough for this case. It creates vector displacement unknowns and element-wise stress/strain storage, then later receives the material model, body force, and boundary-condition handler.

```Python
    # Define momentum equation
    mom_eq = sf.LinearMomentum(grid, theta=0.0)
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
    salt_density = 2169.11
    rho = salt_density*to.ones(mom_eq.n_elems, dtype=to.float64)
    mat.set_density(rho)
```

### Block 13: Constitutive model

The `Spring` element supplies the isotropic elastic stiffness matrix used as the backbone of the mechanical model. Other inelastic mechanisms evolve around this elastic response.

```Python
    # Constitutive model
    E0 = 24.35 * GPa*to.ones(mom_eq.n_elems, dtype=to.float64)
    nu0 = 0.18 * to.ones(mom_eq.n_elems, dtype=to.float64)
    spring_0 = sf.Spring(E0, nu0, "spring")
```

### Block 14: Create creep (A,n)

`DislocationCreep` adds a temperature- and stress-dependent salt-creep mechanism. The parameters `A`, `Q`, and `n` define the Arrhenius-type power law used to update non-elastic strain rates inside the simulator loop.

```Python
    # Create creep (A,n)
    A = 4.51e-33 * to.ones(mom_eq.n_elems, dtype=to.float64)
    Q = 104747.93 * to.ones(mom_eq.n_elems, dtype=to.float64)
    n = 5.5 * to.ones(mom_eq.n_elems, dtype=to.float64)
    creep_0 = sf.DislocationCreep(A, Q, n, "creep_0")
```

### Block 15: Create constitutive model

These calls register elastic stiffness, history-dependent inelastic mechanisms with the `Material` container. Registration order matters because the simulator loops over these lists when assembling tangents and updating internal variables.

```Python
    # Create constitutive model
    mat.add_to_elastic(spring_0)
    if case_name == "creep":
        mat.add_to_non_elastic(creep_0)
```

### Block 16: Set constitutive model

Attaching the material to the momentum equation initializes stiffness-related fields and makes the constitutive mechanisms available to the mechanical solver.

```Python
    # Set constitutive model
    mom_eq.set_material(mat)
```

### Block 17: Set body forces

The gravity vector serves two purposes in mechanics: it builds the volume body-force term through density, and it is passed to pressure boundary objects so hydrostatic variation has the correct vertical sign.

```Python
    # Set body forces
    g = -9.81
    # g = 0.0
    g_vec = [0.0, 0.0, g]
    mom_eq.build_body_force(g_vec)
```

### Block 18: Set initial temperature field

The temperature arrays define the reference and current thermal state used by temperature-dependent mechanical mechanisms.

```Python
    # Set initial temperature field
    T0_field = 326 * to.ones(mom_eq.n_elems)
    mom_eq.set_T0(T0_field)
    mom_eq.set_T(T0_field)
```

### Block 19: Continuation of the script

A parabolic time controller is chosen to distribute a fixed number of time levels nonuniformly between the initial and final times. SafeInCave stores these times internally in seconds, so the examples can be written in hours, days, or years while the solvers use consistent units.

```Python
    tc_eq = sf.TimeControllerParabolic(n_time_steps=20,
										initial_time=0.0,
										final_time=1,
										time_unit="day")
```

### Block 20: Equilibrium Boundary Conditions

These displacement constraints are expressed component by component. SafeInCave locates the DOFs on the named mesh boundary and pins only the requested component, which is how roller or symmetry conditions are represented without overconstraining the model.

```Python
    # -------------------------
    # Equilibrium Boundary Conditions
    # -------------------------
    bc_bottom = momBC.DirichletBC(
        boundary_name="BOTTOM",
        component=2,
        values=[0.0, 0.0],
        time_values=[0.0, tc_eq.t_final]
    )
    bc_west = momBC.DirichletBC(
        boundary_name="WEST",
        component=0,
        values=[0.0, 0.0],
        time_values=[0.0, tc_eq.t_final]
    )
    bc_east = momBC.DirichletBC(
        boundary_name="EAST",
        component=0,
        values=[0.0, 0.0],
        time_values=[0.0, tc_eq.t_final]
    )
    bc_south = momBC.DirichletBC(
        boundary_name="south",
        component=1,
        values=[0.0, 0.0],
        time_values=[0.0, tc_eq.t_final]
    )
    bc_north = momBC.DirichletBC(
        boundary_name="north",
        component=1,
        values=[0.0, 0.0],
        time_values=[0.0, tc_eq.t_final]
    )
```

### Block 21: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
    bc_equilibrium = momBC.BcHandler(mom_eq)
    bc_equilibrium.add_boundary_condition(bc_bottom)
    bc_equilibrium.add_boundary_condition(bc_west)
    bc_equilibrium.add_boundary_condition(bc_east)
    bc_equilibrium.add_boundary_condition(bc_south)
    bc_equilibrium.add_boundary_condition(bc_north)
```

### Block 22: Continuation of the script

The boundary handler is the bridge between simple Python boundary objects and the assembled variational problem. It stores each condition by type, updates time-dependent values, and supplies the forms or DOF constraints needed by the solver.

```Python
    mom_eq.set_boundary_conditions(bc_equilibrium)
```

### Block 23: Output (equilibrium)

This path groups all files from the run. `SaveFields` later creates one XDMF subfolder per registered field, and cavern handlers use the same folder to store fluid-operation histories.

```Python
    # Output (equilibrium)
    output_mom = sf.SaveFields(mom_eq)
    output_mom.set_output_folder(os.path.join(output_folder, "equilibrium"))
    output_mom.add_output_field("u", "Displacement (m)")
    output_mom.add_output_field("eps_tot", "Total strain (-)")
    output_mom.add_output_field("sig", "Stress (Pa)")
    output_mom.add_output_field("p_elems", "Mean stress (Pa)")
    output_mom.add_output_field("q_elems", "Von Mises stress (Pa)")
    outputs = [output_mom]
```

### Block 24: Continuation of the script

`Simulator_M` advances the mechanical problem only. When a cavern handler is supplied, the simulator still updates cavern volumes and pressure boundary loads, but it does not solve a surrounding heat equation.

```Python
    sim = sf.Simulator_M(mom_eq, tc_eq, outputs, compute_elastic_response=True)
    sim.run()
```

### Block 25: run("creep")

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python

def main():
    run("elastic")
    # run("creep")
```

### Block 26: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python
if __name__ == '__main__':
    main()
```

## Main Functions

- `safeincave.GridHandlerGMSH`: loads the settlement mesh from the example grid folder.
- `safeincave.LinearMomentum`: defines the mechanical balance equation.
- `safeincave.Material`: stores density and constitutive model data.
- `safeincave.Spring`: defines the elastic stiffness.
- `safeincave.DislocationCreep`: adds time-dependent creep deformation.
- `safeincave.TimeControllerParabolic`: creates the long-term time discretization.
- `safeincave.MomentumBC.DirichletBC`: constrains the model boundaries.
- `safeincave.SaveFields`: saves displacement and stress-related fields.
- `safeincave.Simulator_M`: runs the settlement-oriented mechanical simulation.

## Running the Example

From the repository root, run:

```bash
cd examples/mechanics/8_settlement
python main.py
```

The folder also contains `main_settlement_0.py`, which can be used for an alternative settlement setup:

```bash
python main_settlement_0.py
```

Results are written to `output/`.

## Conclusion

This example demonstrates how SafeInCave can be used for settlement-oriented mechanical calculations. The saved displacement fields can be post-processed to quantify surface movement above salt caverns or other subsurface features.
