# Example 3: Mass Exchange

## Introduction

This example computes the thermodynamic response of hydrogen during injection and withdrawal in a constant-volume cavern. The imposed mass increments drive changes in pressure, temperature, and density.

## Problem description

This section walks through `examples/thermodynamics/3_mass_exchange/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

### Block 1: Imports and dependencies

The example stays at the lumped-fluid level, so it imports numerical arrays, plotting tools, and `CavernThermodynamics` rather than finite-element classes. `CavernThermodynamics` is the SafeInCave wrapper around CoolProp and the nonlinear energy/mass balance used to advance the fluid state.

```Python
import numpy as np
import matplotlib.pyplot as plt
from safeincave.Thermodynamics import CavernThermodynamics
```

### Block 2: Define conversion units

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
# Define conversion units
minute = 60
hour = 60*minute
day = 24*hour
cm = 1e-2
MPa = 1e6
ton = 1000
Mton = ton*1e6
```

### Block 3: Function `apply_grey_theme`

This is a plotting convenience, not a modeling operation. It gives all thermodynamic plots a common visual style so differences in pressure, temperature, or density are easier to compare.

```Python
def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color='0.85'):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color=grid_color)
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor(back_color)
```

### Block 4: Function `mass_rate`

This helper creates the injection or withdrawal history. Positive increments add hydrogen to the fixed-volume cavern, while negative values remove mass.

```Python


def mass_rate(M_min, M_max):
    time = np.linspace(0, 100, 100)
    dM_list = M_min + (M_max - M_min)*time/time[-1]
    return time, dM_list
```

### Block 5: Function `calculate_thermodynamics`

This function contains the physical update loop. It creates a `CavernThermodynamics` object for the requested fluid, computes the initial density/internal-energy state, and repeatedly calls `solve` so mass exchange is translated into pressure, temperature, and density.

```Python


def calculate_thermodynamics(fluid_name, dM_list, P_init=10*MPa, T_init=300):
    model = CavernThermodynamics(fluid_name)
    P_atm = 101325.0
    rho_init, _, _ = model.rho_u_h(P_init + P_atm, T_init)
    T, P, rho = [T_init], [P_init], [rho_init]
    for i in range(1, len(dM_list)):
        P1, T1, rho1 = model.solve( dm = dM_list[i],
                                    Q_in = 0.0,
                                    T_in = 280.0,
                                    P0 = P[i-1] + P_atm,
                                    T0 = T[i-1],
                                    V0 = 1,
                                    V1 = 1)
        P.append(P1 - P_atm)
        T.append(T1)
        rho.append(rho1)
    return np.array(P), np.array(T), np.array(rho)
```

### Block 6: Function `main`

This function is the procedural spine of the example. The following blocks fill in the mesh, equation objects, material model, boundary data, outputs, and simulator call that this function orchestrates.

```Python


def main():
    fig, axis = plt.subplots(2, 4, figsize=(16, 8))
    fig.subplots_adjust(top=0.935, bottom=0.155, left=0.062, right=0.980, hspace=0.35, wspace=0.38)
    ax1, ax2, ax3, ax4 = axis[0,:]
    ax5, ax6, ax7, ax8 = axis[1,:]
```

### Block 7: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
    time, dM_in = mass_rate(0.5, 0.5)
    P, T, rho = calculate_thermodynamics("Hydrogen", dM_in)
```

### Block 8: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax1.plot(time, dM_in, "-", linewidth=2.0, color="black")
    ax1.set_xlabel("Time (s)", fontname="serif", size=12)
    ax1.set_ylabel("Mass (kg)", fontname="serif", size=12)
```

### Block 9: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax2.plot(time, P/MPa, "-", linewidth=2.0, color="steelblue", label=r"H$_2$")
    ax2.set_xlabel("Time", fontname="serif", size=12)
    ax2.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
```

### Block 10: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax3.plot(time, T-273, "-", linewidth=2.0, color="steelblue", label=r"H$_2$")
    ax3.set_xlabel("Time", fontname="serif", size=12)
    ax3.set_ylabel("Temperature (°C)", fontname="serif", size=12)
```

### Block 11: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax4.plot(time, rho, "-", linewidth=2.0, color="steelblue", label=r"H$_2$")
    ax4.set_xlabel("Time", fontname="serif", size=12)
    ax4.set_ylabel("Density (kg/m3)", fontname="serif", size=12)
```

### Block 12: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

    time, dM_out = mass_rate(-0.05, -0.05)
    P, T, rho = calculate_thermodynamics("Hydrogen", dM_out)
```

### Block 13: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax5.plot(time, dM_out, "-", linewidth=2.0, color="black")
    ax5.set_xlabel("Time (s)", fontname="serif", size=12)
    ax5.set_ylabel("Mass (kg)", fontname="serif", size=12)
```

### Block 14: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax6.plot(time, P/MPa, "-", linewidth=2.0, color="steelblue", label=r"H$_2$")
    ax6.set_xlabel("Time", fontname="serif", size=12)
    ax6.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
```

### Block 15: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax7.plot(time, T-273, "-", linewidth=2.0, color="steelblue", label=r"H$_2$")
    ax7.set_xlabel("Time", fontname="serif", size=12)
    ax7.set_ylabel("Temperature (°C)", fontname="serif", size=12)
```

### Block 16: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax8.plot(time, rho, "-", linewidth=2.0, color="steelblue", label=r"H$_2$")
    ax8.set_xlabel("Time", fontname="serif", size=12)
    ax8.set_ylabel("Density (kg/m3)", fontname="serif", size=12)
```

### Block 17: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    apply_grey_theme(fig, axis.flatten())
```

### Block 18: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    plt.show()
```

### Block 19: Script entry point

The entry-point guard keeps the file importable while still running `main()` when executed as a script.

```Python


if __name__ == "__main__":
	main()
```

## Main Functions

- `safeincave.Thermodynamics.CavernThermodynamics`: solves the cavern-fluid thermodynamic balance.
- `CavernThermodynamics.rho_u_h`: initializes the hydrogen state.
- `CavernThermodynamics.solve`: updates the state using the prescribed mass exchange `dm` and inlet temperature `T_in`.
- `mass_rate`: creates injection and withdrawal histories.
- `calculate_thermodynamics`: iterates through the mass history and records pressure, temperature, and density.

## Running the Example

From the repository root, run:

```bash
cd examples/thermodynamics/3_mass_exchange
python main.py
```

The script opens a Matplotlib figure showing mass exchange and the resulting hydrogen state variables.

## Conclusion

This example demonstrates how `CavernThermodynamics` handles mass exchange. The workflow is directly applicable to gas-storage operation schedules where injection and withdrawal rates are prescribed.
