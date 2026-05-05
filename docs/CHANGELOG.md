# Changelog

## Next release
- Implemented stabilized mixed formulation. Unknowns are displacement and mean stress fields
- Included proper tests for linear elasticity model
- Fixed issue of Robin boundary condition not being updated
- Implemented new classes for handling caverns with different operational conditions
- Implemented thermodynamic model for fluid (gas/liquid) inside cavern
- Implemented initial stresses to both P1 and P1-P1 formulations

## 2.0.0
- Implemented MPI parallelisation
- Implemented heat diffusion equation for thermal effects, including thermal strains
- Implemented pressure solution creep model
- Changed output format to XDMF
