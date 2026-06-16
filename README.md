<img src="https://raw.githubusercontent.com/ADMIRE-Public/SafeInCave/develop/assets/logo_2.png" width="500"/>

[![Version](https://img.shields.io/badge/version-3.0.3-blue)](https://github.com/ADMIRE-Public/SafeInCave)
[![Platform](https://img.shields.io/badge/Platform-Ubuntu%20%7C%20Windows%20(WSL)-blue)](https://ubuntu.com/wsl)  
[![FEniCSx](https://img.shields.io/badge/Dependency-FEniCSx%200.9.0-important)](https://fenicsproject.org)


## Overview
SafeInCave is an open-source 3D finite element simulator specifically designed for salt cavern mechanics. It combines the most relevant physical phenomena with a robust numerical scheme, placing SafeInCave as a powerful tool for performance assessments of salt caverns under different operational conditions. 

<img src="https://raw.githubusercontent.com/ADMIRE-Public/SafeInCave/develop/docs/images/multicavern_q.gif" width="500"/>


## Open-source
Being **open-source** is one of the most important aspect of SafeInCave. This has important implications for both academic and industrial applications. Full transparency, auditability, and reproducibility are essential for building trust in simulation results. Moreover, access to the source code allows users to inspect, adapt, and extend the simulator when needed, enabling contributions that can benefit the entire salt community.

## Key Features

- **Tetrahedral meshes**: Easily discretize complex geometries without compromising numerical stability.
- **Thermodynamics**: Different fluids with different operational conditions are possible.
- **Coupled physics**: Mechanics, heat diffusion, and thermodynamic models are fully coupled.
- **Thermal effects**: Solve heat diffusion equation and include thermal strains and creep thermal responses.
- **Constitutive model**: Munson-Dawson model, two branches model, thermal strains, Cam-Clay model, etc.
- **Robust linearization**: Provides robustness and flexibility to include new constitutive models.
- **Time discretization**: Choose between Explicit, Crank-Nicolson, and Fully-Implicit schemes.
- **XDMF output**: Efficient output format in terms of size and postprocessing.


## Installation
SafeInCave installation depends on [FEniCSx](https://fenicsproject.org/) installation. For Windows users, the installaion pipeline consists of:

1) Installing [WSL](https://learn.microsoft.com/en-us/windows/wsl/)

2) Installing Ubuntu

3) Installing [FEniCSx](https://fenicsproject.org/download/)

4) Installing [SafeInCave](https://admire-public.github.io/SafeInCave/home/installation/)

See SafeInCave [documentation](https://admire-public.github.io/SafeInCave/) for full documentation (under active development).


## Getting started
Users can build their own simulators using the *safeincave* package. The [documentation](https://admire-public.github.io/SafeInCave/) webpage shows detailed examples of how to set up purely **mechanical**, **heat diffusion**, and **thermomechanical** simulations. These examples show how to build constitutive models, apply different types of boundary conditions, assign material properties, etc.


## Mantainers
- [Davi R. Damasceno](https://www.linkedin.com/in/drdamasceno/)
- [Gijs van den Brekel](https://www.linkedin.com/in/gijs-van-den-brekel-041866229/)
- [Hadi Hajibeygi](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-engineering/sections/reservoir-engineering/staff/academic-staff/profdr-h-hajibeygi)
- [Hermínio Tasinafo Honório](https://www.linkedin.com/in/herminioth/)
- [Lucas Landeweerd](https://www.linkedin.com/in/lucaslandeweerd/)
- [Mathias Jose Kreutz Erdtmann](https://www.linkedin.com/in/mathias-kreutz-erdtmann/)


## License
SafeInCave is released under the **BSD 3-Clause License** (BSD-3-Clause), a permissive open-source license that allows use, modification, and redistribution in both open-source and proprietary software.

See the [LICENSE](https://opensource.org/license/BSD-3-clause) file for the complete license terms.


## Key publications
- Honório, H.T, Hajibeygi, H. [Three-dimensional multi-physics simulation and sensitivity analysis of cyclic hydrogen storage in salt caverns](https://doi.org/10.1016/j.ijhydene.2024.11.081). International Journal of Hydrogen Energy, 2024.

- Honório, H.T., Franceschini, A., Ferronato, M., & Hajibeygi, H. [Salt cavern simulations with a stabilized mixed finite element formulation for low-order tetrahedral elements](https://doi.org/10.1016/j.cma.2026.119073). CMAME, 2026.

- Honório, H.T., Amini, M.S., Landeweerd, L., & Hajibeygi, H. [SafeInCave: An Open-Source Simulator for Energy Storage in Heterogeneous Salt Caverns](https://smri.memberclicks.net/assets/docs/Fall2025/TechSessions/20_MP2025F_Honorio.pdf). In Proceedings of the SMRI Fall Technical Conference, Whichita, Kansas, US, 2025.


## Acknowledgements
- [Shell Global Solutions International B.V](https://www.shell.com/) for sponsoring the [project SafeInCave](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-engineering/research/research-themes/energy-transition/safeincave), within which this simulator was developed.
- [Energi Simulation](https://energisimulation.com/) for currently supporting this project.

