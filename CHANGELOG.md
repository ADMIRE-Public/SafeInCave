# Changelog

## 2.1.0
- Implemented stabilized mixed formulation. Unknowns are displacement and mean stress fields
- Included proper tests for linear elasticity model
- Fixed issue of Robin boundary condition not being updated

## 2.0.0
- Implemented MPI parallelisation
- Implemented heat diffusion equation for thermal effects, including thermal strains
- Implemented pressure solution creep model
- Changed output format to XDMF
