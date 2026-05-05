# Example 2: Isochoric Heat Exchange

## Introduction

This example studies heat addition and heat removal in a constant-volume cavern. The fluid is water, and the model computes how pressure, temperature, and density evolve when heat is exchanged without mass transfer or volume change.

## Problem description

This section walks through `examples/thermodynamics/2_isochoric_heat_exchange/main.py` as an annotated, notebook-style listing. The code blocks follow the script in order and include all non-empty lines from `main.py`. Each block explains not only what the code does, but why the SafeInCave object is needed in the simulation workflow.

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

### Block 4: Function `heat`

This helper creates the heat-input history. The subsequent thermodynamic loop keeps volume and mass fixed, so `Q_in` is the only source term changing the water state.

```Python


def heat(Q_min, Q_max):
    time = np.linspace(0, 100, 100)
    Q_list = Q_min + (Q_max - Q_min)*time/time[-1]
    return time, Q_list
```

### Block 5: Function `calculate_thermodynamics`

This function contains the physical update loop. It creates a `CavernThermodynamics` object for the requested fluid, computes the initial density/internal-energy state, and repeatedly calls `solve` so heat input is translated into pressure, temperature, and density.

```Python


def calculate_thermodynamics(fluid_name, Q_list, P_init=10*MPa, T_init=300):
    model = CavernThermodynamics(fluid_name)
    P_atm = 101325.0
    rho_init, _, _ = model.rho_u_h(P_init + P_atm, T_init)
    T, P, rho = [T_init], [P_init], [rho_init]
    for i in range(1, len(Q_list)):
        P1, T1, rho1 = model.solve( dm = 0.0,
                                    Q_in = Q_list[i],
                                    T_in = 0.0,
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
    fig.subplots_adjust(top=0.935, bottom=0.155, left=0.05, right=0.980, hspace=0.35, wspace=0.38)
    ax1, ax2, ax3, ax4 = axis[0,:]
    ax5, ax6, ax7, ax8 = axis[1,:]
```

### Block 7: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python
    time, Q_heat_up = heat(7e5, 7e5)
    P, T, rho = calculate_thermodynamics("Water", Q_heat_up)
```

### Block 8: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax1.plot(time, Q_heat_up/1e3, "-", linewidth=2.0, color="black")
    ax1.set_xlabel("Time (s)", fontname="serif", size=12)
    ax1.set_ylabel("Heat (kJ)", fontname="serif", size=12)
```

### Block 9: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax2.plot(time, P/MPa, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax2.set_xlabel("Time", fontname="serif", size=12)
    ax2.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
```

### Block 10: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax3.plot(time, T-273, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax3.set_xlabel("Time", fontname="serif", size=12)
    ax3.set_ylabel("Temperature (°C)", fontname="serif", size=12)
```

### Block 11: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax4.plot(time, rho, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax4.set_xlabel("Time", fontname="serif", size=12)
    ax4.set_ylabel("Density (kg/m3)", fontname="serif", size=12)
    ax4.set_ylim(1000.4, 1001.4)
```

### Block 12: Continuation of the script

This small block connects the surrounding setup steps and is retained so the annotated listing remains complete and in the same order as `main.py`.

```Python

    time, Q_cool_down = heat(-7e5, -7e5)
    P, T, rho = calculate_thermodynamics("Water", Q_cool_down)
```

### Block 13: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax5.plot(time, Q_cool_down/1e3, "-", linewidth=2.0, color="black")
    ax5.set_xlabel("Time (s)", fontname="serif", size=12)
    ax5.set_ylabel("Heat (kJ)", fontname="serif", size=12)
```

### Block 14: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax6.plot(time, P/MPa, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax6.set_xlabel("Time", fontname="serif", size=12)
    ax6.set_ylabel("Pressure (MPa)", fontname="serif", size=12)
```

### Block 15: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax7.plot(time, T-273, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax7.set_xlabel("Time", fontname="serif", size=12)
    ax7.set_ylabel("Temperature (°C)", fontname="serif", size=12)
```

### Block 16: Continuation of the script

This plotting block presents the prescribed input and the resulting state variables. The thermodynamics examples use plots as their primary output rather than writing XDMF fields.

```Python
    ax8.plot(time, rho, "-", linewidth=2.0, color="lightcoral", label=r"H$_2$O")
    ax8.set_xlabel("Time", fontname="serif", size=12)
    ax8.set_ylabel("Density (kg/m3)", fontname="serif", size=12)
    ax8.set_ylim(1000.4, 1001.4)
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

- `safeincave.Thermodynamics.CavernThermodynamics`: provides fluid property evaluation and thermodynamic state updates.
- `CavernThermodynamics.rho_u_h`: initializes the water state from pressure and temperature.
- `CavernThermodynamics.solve`: advances the cavern state using `Q_in` while keeping `V0` and `V1` equal.
- `heat`: creates the heat-input histories.
- `calculate_thermodynamics`: applies the heat history and collects the resulting state variables.

## Running the Example

From the repository root, run:

```bash
cd examples/thermodynamics/2_isochoric_heat_exchange
python main.py
```

The script displays Matplotlib plots for the heat history and the resulting pressure, temperature, and density response.

## Conclusion

This example isolates the thermal part of cavern-fluid thermodynamics. It is useful for checking heat-transfer assumptions before combining the fluid model with heat diffusion through the rock.
