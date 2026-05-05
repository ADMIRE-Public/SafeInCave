# Thermodynamic model
The thermodynamic models implemented in SafeInCave use [CoolProp](https://coolprop.org/) for the equations of state. Substances of interest, such as methane, hydrogen, air, and water, are supported. Brine is treated in CoolProp as an incompressible fluid, it is currently not supported.  

In SafeInCave, the operational conditions of caverns can be specified in two ways:

- Specified pressure ($P$) and temperature ($T$) histories;
- Specified mass flow rate ($\dot{m}$) history.

These two cases require different thermodynamic model implementations and they are described below.

## Imposed pressure-temperature history
Specifying pressure and temperature over time means that the thermodynamic state of the substance is already fully pre-defined. Therefore, for this case, the only responsability of the thermodynamic model is to calculate the fluid's density over time. This calculation is performed automatically in SafeInCave.

Using CoolProp, density can be calculated based on P and T as shown in the snippet below.

```python
import CoolProp.CoolProp as CP
# Specify fluid, P, and T
AS = CP.AbstractState("HEOS", fluid)
AS.update(CP.PT_INPUTS, P, T)
rho = self.AS.rhomass()
```

Although it is considered perfect mixing of the fluid such that temperature is homogeneous, pressure still varies with depth according to fluid's specific weight, hence the importance of calculating density for each $P-T$ pair.

>**Note:** Density is considered to be constant inside the cavern, even though it should follow the vertical pressure variation.

In addition to imposing temperature $T^*$ and pressure $P^*$, the user also needs to inform the position $z^*$ where this pressure refers to. In this manner, the vertical pressure variation is calculated as

$$
p(z) = P^* + \rho g (z^* - z)
$$

with the $z$ axis pointing downwards, and $\rho=\rho(P^*, T^*)$.


## Imposed mass flow rate history
In this case, the user spcifies the initial temperature ($T_i$) and pressure $P_i$ (hence, thermodynamic state is defined), and the mass flow rate history. Therfore, pressure, temperature, and density histories are a result of the thermodynamic model.  

Note that the initial density is defined as $\rho_i = \rho(P_i, T_i)$. Additionally, the initial cavern volume $V_i$ is also known, and the initial mass is determined as $M_i = \rho_i V_i$. Therefore, with $M_i$ and the history of mass flow rates fully determines the mass history inside the cavern.

The $P-T$ histories are obtained from the first law of thermodynamics. Considering the cavern as an open system, the first law of thermodynamics states that variation of internal energy between states 0 and 1 is given by

$$
U_1 - U_0 = Q - W + M_\text{in} h_\text{in} - M_\text{out} h_\text{out}
\label{eq:first_law_0}
$$

where $U_1 - U_0$ represent the internal energy variation between states (or time steps) 0 and 1 , $Q$ are the heat transferred through the cavern walls between states 0 and 1, $W$ is the work done by the cavern walls between 0 and 1, $M_\text{in/out}$ is the total mass that injected/removed from the cavern between 0 and 1, and $h_\text{in/out}$ denotes the enthalpy of the fluid and entered/exited the cavern. This is illustrated in [](#fig-cavern-thermo)-a and [](#fig-cavern-thermo)-b.

![(a) Cavern withdrawal, (b) cavern injection, and (c) heat transfer through cavern walls.](../../images/caverns.svg){#fig-cavern-thermo width="110%"}


### Work calculation
The work done by the cavern walls on the system is given by

$$
W = \int_{V_0}^{V_1} p \mathrm{d}V = \bar{p} \Delta V,
$$

where $\Delta V$ is the cavern volume variation between states 0 and 1. Additionally, $\bar{p}$ is the pressure acting on the cavern walls between the two states, which is considered to be

$$
\bar{p} = P_1.
\label{eq:p_bar_W}
$$

>**Note 1:** The pressure assumption in Eq. $\eqref{eq:p_bar_W}$ is equivalent to a fully-implicit time integration approach.

>**Note 2:** From the thermodynamic's model perspective, $V_1$ is a given, but it's actually calculated by the mechanical model in SafeInCave.


### Total heat transfer
The total heat transferred through the cavern walls is calculated based on the temperature gradients obtained from the heat diffusion model, as illustrated in [](#fig-cavern-thermo)-c. That is,

$$
Q = \int_{t_0}^{t_1} \dot{Q} \mathrm{d}t = \dot{Q}_1 \Delta t, 
\quad \text{where} \quad
\dot{Q}_1 = \oint_\Gamma -k \nabla T \cdot \hat{\mathbf{n}} \mathrm{d}\Gamma
\label{eq:Q_0}
$$

where $\dot{Q}_1$ is the heat transfer rate (in $W$) calculated with the gradient of the most updated temperature field $\nabla T$ and the heat conduction coefficient $k$. Moreover, the closed integral in Eq. $\eqref{eq:Q_0}$ represent is performed on the entire cavern wall $\Gamma$, which outward oriented by the unitary normal vector $\hat{\mathbf{n}}$ (see [](#fig-cavern-thermo)-c).

### Energy balance
In general, the energy balance equation for the salt cavern can be expressed as

$$
M_1 u(P_1, T_1) - M_0 u(P_0, T_0) = Q - W + M_\text{in/out} h(P_\text{in/out}, T_\text{in/out}).
\label{eq:energy_balance_0}
$$

where $u$ is the specific internal energy and $M_\text{in/out} = M_1 - M_0$. Note that:

- If $M_1 > M_0$, then $M_\text{in/out} = M_\text{in} > 0$ (injection period);
- If $M_1 < M_0$, then $M_\text{in/out} = M_\text{out} < 0$ (production period);
- If $M_1 = M_0$, then $M_\text{in/out} = 0$ (stand still period);

During injection period, the user has to specify temperature $T_\text{in}$ of the injected fluid. In order to define the thermodynamic state of the injected fluid, the pressure $P_\text{in}$ must also be known. Due to the difficulty in specifying this pressure, it is calculate internally as

$$
P_\text{in} = \frac{P_0 + P_1}{2}.
$$

During production period, the temperature and pressure of the fluid exiting the cavern are computed as

$$
T_\text{out} = \frac{T_0 + T_1}{2} \quad \text{and} \quad P_\text{out} = \frac{P_0 + P_1}{2}.
$$

### Mass balance

Note that the thermodynamic's model perspective, the total heat tranfer $Q$, the work $W$, and the final volume $V_1$ are inputs, just like $P_0$, $T_0$, $M_0$, and $M_1$. The only two unknowns are $P_1$ and $T_1$. The closure equation is the mass balance, which states that

$$
\rho(P_1, T_1) = \frac{M_1}{V_1}.
\label{eq:mass_balance_0}
$$


### Final model equations
The solution of the thermodynamic model consists of finding $P_1$ and $T_1$ such that both Equations $\eqref{eq:energy_balance_0}$ and $\eqref{eq:mass_balance_0}$ are satisfied, as illustrated in [](#fig-thermodynamic-model). This a non-linear problem to solve since $u=u(P,T)$, $h=h(P,T)$, and $\rho=\rho(P,T)$ are all non-linear functions. A Newton-like approach is employed to solve this non-linear system of equations. 

![Thermodynamic model: inputs, model equations and outputs.](../../images/thermodynamic_box.svg){#fig-cavern-thermo width="90%"}




