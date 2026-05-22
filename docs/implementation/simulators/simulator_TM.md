# Thermo-mechanical simulator (TM)
This simulator is designed to handle thermo-mechanical problems, in which the cavern mechanics is also affected by temperature distribution. Consequently, it requires objects of both momentum equation (either *LinearMomentum* or *LinearMomentumMixed*) and heat diffusion (*HeatDiffusion*). As the other simulators, *Simulator_TM* also receives the usual *TimeController*, *Outputs*, and *CavernHandler* objects.

???+ example "Example"
    ```python
    sim = Simulator_TM(eq_mom, eq_heat, time_control, outputs, caverns)
    ```

The *CavernHandler* object is optional. If a *None* object is provided instead, then the user needs to provide appropriate boundary conditions for the cavern walls in both mechanical and heat diffusion models. If *CavernHandler* is provided, then the cavern walls are subjected to a dynamic *Neumann* boundary condition for the momentum equation, and *Robin* boundary condition is imposed on the heat diffusion equation.

The thermo-mechanical coupling in salt cavern is largely adopted in the literature. The fact that salt deformations occur very slowly, the heat generate by such deformations is fairly negligible. As a result, the mechanical model has no effect on the heat diffusion process. On the other side, the temperature evolution has a strong impact on the mechanical behavior of the salt caverns as it (i) causes thermal strains to appear and (ii) strongly affect dislocation and pressure solution creep through the Arrhenius law.

The previous paragraph suggests a **one-way** coupling between mechanical and thermal models. However, the presence of a thermodynamic model for the fluid inside the cavern introduces a **two-way** coupling between the two governing equations. This is because the temperature and pressure of the fluid inside the cavern -- used as boundary conditions for mechanical and thermal model, respectively -- are a result of the thermodynamic model, which depends on (i) the heat transferred through the cavern walls, and (ii) the work done by the cavern walls. The heat transferred through the cavern walls to the fluid inside the cavern is calculated from heat diffusion model. The work done by the cavern walls on the cavern fluid is calculated by mechanical model. With heat $Q$ and work $W$, the thermodynamic model is calculated and provides temperature ($T_1$), pressure ($P_1$), and density ($\rho_1$) for the fluid. These quantities are used as boundary conditions for the two governing equations, as illustrated in [](#fig-simulator-TM). 

This process is repeated until convergence is reached, before advancing to the next time level. Notice that these iterations simultaneously handle the coupling between the governing equations and the non-linearities of the mechanical model.

![Coupling scheme for the mechanical simulator.](../../images/simulator_TM.svg){#fig-simulator-TM .wide-figure style="width: min(700px, 95vw);"}