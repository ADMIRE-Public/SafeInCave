# Changelog

## Next release
- Implemented stabilized mixed formulation. Unknowns are displacement and mean stress fields
- Included proper tests for linear elasticity model
- Fixed issue of Robin boundary condition not being updated
- Implemented new classes for handling caverns with different operational conditions
- Implemented thermodynamic model for fluid (gas/liquid) inside cavern
- Implemented initial stresses to both P1 and P1-P1 formulations
- Included firs draft for documentation
- Implemented CAM-CLAY model
- Implemented Munson-Dawson model
- Fixed bug in Desai's model related to sigma_t
- Removed dependency of BcHandlers from Equation
- Deleted SimulatorFull class
- Fixed #19: Standardized compute_CT input arguments across MomentumEquation classes.
- Implemented solver definition internally in the Equations
- ML model for h parameter calculation now uses skops to be loaded.
- Changed license to BSD-3-Clause.

## 2.0.0
- Implemented MPI parallelisation
- Implemented heat diffusion equation for thermal effects, including thermal strains
- Implemented pressure solution creep model
- Changed output format to XDMF
