<!-- # SafeInCave Documentation -->
# 

![](images/logo_2.png){#fig-constitutive-model-0 width="80%"}

[![](https://img.shields.io/badge/version-2.1.0-blue)](https://github.com/ADMIRE-Public/SafeInCave)
[![](https://img.shields.io/badge/DOI-10.1234%2Fabcde.12345678-blue)](https://doi.org/10.5281/zenodo.17169591)
[![](https://img.shields.io/badge/Platform-Ubuntu%20%7C%20Windows%20(WSL)-blue)](https://ubuntu.com/wsl)  
[![](https://img.shields.io/badge/Dependency-FEniCSx%200.9.0-important)](https://fenicsproject.org)
[![](https://img.shields.io/github/license/ADMIRE-Public/SafeInCave)](https://www.gnu.org/licenses/gpl-3.0.en.html)



<!-- Welcome to the SafeInCave documentation page. SafeInCave is an open-source simulator for the thermo-mechanical behavior of salt caverns.  -->

<!-- --- -->

## Overview
SafeInCave is an open-source 3D finite element simulator specifically designed for salt caverns. It combines the most relevant physical phenomena with a robust numerical scheme, making SafeInCave as a powerful tool for performance assessments of salt caverns under different operational conditions. 

<!-- The thermodynamic model for the cavern fluid is fully integrated with the mechanical and thermal behavior of the cavern. In this manner, it can seemlessly simulate both storage and abandonment phases considering different types of fluids (water, hydrogen, methane, etc). SafeInCave also features advanced constitutive models for salt rocks, such as the Munson-Dawson model. Finally, it is the only simulator able to provide numerically stable simulations with tetrahedral meshes, which ensure great geometrical flexibility. -->

## Open-source
Being **open-source** is perhaps the most important aspect of SafeInCave. This has many implications for both academic and industrial applications. Being fully transparent, auditable and reproducible is essencial to build trust in the results. Moreover, having access to the source code and being able to modify it as necessary allows for anyone to make contributions that can benefit the whole salt community!

<!-- --- -->

## Key Features

<!-- - **MPI-powered parallelism**: Scale simulations efficiently with mpi4py for distributed computing -->
- **Tetrahedral meshes**: Easily discretize complex geometries without compromising numerical stability
- **Thermodynamics**: Different fluids with different operational conditions are possible
- **Coupled physics**: Mechanics, heat diffusion, and thermodynamic models are fully coupled
- **Thermal effects**: Solve heat diffusion equation and include thermal strains and creep thermal responses
- **Constitutive model**: Munson-Dawson model, two branches model, thermal strains, Cam-Clay model
- **Robust linearization**: Provides robustness and flexibility to include new constitutive models
- **Time discretization**: Choose between Explicit, Crank-Nicolson, and Fully-Implicit schemes
- **XDMF output**: Efficient output format in terms of size and postprocessing

<!-- --- -->

## Extra material
Video lectures and video tutorials can be found in the [ADMIRE](https://www.youtube.com/@ADMIRE1/featured) YouTube channel. The following videos are currently available:

1) [Tensorial operations (theory)](https://youtu.be/w5KX3F_rdzU?si=QQLVBq1NcrvOiS32)

2) [Tensorial operations (exercises)](https://www.youtube.com/watch?v=JiN6jwp0RPk&t=0s)

3) [Constitutive modeling](https://www.youtube.com/watch?v=fCeJIbjIL10)

4) Stay tuned to [ADMIRE](https://www.youtube.com/@ADMIRE1/featured) YouTube channel for upcoming video lectures.

<!-- --- -->

## Mantainers
- [Hermínio Tasinafo Honório](https://www.linkedin.com/in/herminioth/)
- [Hadi Hajibeygi](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-engineering/sections/reservoir-engineering/staff/academic-staff/profdr-h-hajibeygi)
- [Mathias Erdtmann](https://www.linkedin.com/in/mathias-kreutz-erdtmann/)

---

## How to cite
If you use SafeInCave in your publications, you can cite it as follows:

ADMIRE-Public. “ADMIRE-Public/SafeInCave: SafeInCave V2.0.0”. Zenodo, September 21, 2025. https://doi.org/10.5281/zenodo.17169591.

---

## Acknowledgements
We would like to thank:
- [Shell Global Solutions International B.V](https://www.shell.com/) for sponsoring the [project SafeInCave](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-engineering/research/research-themes/energy-transition/safeincave), within which this simulator was developed.
- [Energi Simulation](https://energisimulation.com/) for currently supporting this project.