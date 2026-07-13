# Changelog

## Next release
- Refactored MaterialProps.py by moving classes into individual files in ConstitutiveModels folder.
- Import automation of files in ConstitutiveModels so that the only step to add a new constitutive model is to include the file in the folder.
- Fixed relative paths inside pytest so that it is more robust to run from anywhere in the repo.
- Fixed import standards to comply with linting.
- Included closing of output.
- Implemented unified output file (solution.xdmf) containing all simulation fields.
- Fix HDF5 file contention (solution.xdmf) in multi-physics simulations.
- Encapsulated convergence criteria into a separate file.
- Fixed dt display during simulations.
- Implemented results logging at a specific point.
- Fixed type hints in constitutive models.
- Included initial state to results tracking.
- Added a plastic yield-consistency convergence gate.
- Made force-residual criterion Dirichlet-aware and displacement-increment criterion robust to zero-increment steps.
- Option for smoothed output functionality.
- Fixed broken von Mises extractor in SimulationLogging (raised AttributeError) and reduced logging setup verbosity.
- Implemented local extensions mechanism.
- Added a closing solve per time step so the committed inelastic strain uses the final converged rate (fixes hardening drift on models with persistent internal variables); constitutive models can now opt in to an exact elastic-trial reconstruction via uses_exact_trial.

## 3.0.3
- Fixed mean stress calculation in P1P1 (mixed) formulation
- Removed mean stress smoother from P1 formulation

## 3.0.2
- Fixed bug in initial stress implementation in the mixed formulation

## 3.0.1
- Fixed License typos.
- Updated patch version

## 3.0.0
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
- Added automatic conda environment installations.

## 2.0.0
- Implemented MPI parallelisation
- Implemented heat diffusion equation for thermal effects, including thermal strains
- Implemented pressure solution creep model
- Changed output format to XDMF