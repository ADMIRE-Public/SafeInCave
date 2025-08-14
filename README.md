![alt text](<./docs_archive/source/_static/logo_2.png>){width=35%}

[![Version](https://img.shields.io/badge/version-2.0.0-blue)](https://gitlab.tudelft.nl/ADMIRE_Public/safeincave)
[![Platform](https://img.shields.io/badge/Platform-Ubuntu%20%7C%20Windows%20(WSL)-blue)](https://ubuntu.com/wsl)  
[![FEniCSx](https://img.shields.io/badge/Dependency-FEniCSx%200.9.0-important)](https://fenicsproject.org)


> **Note**: This project requires **FEniCSx** and runs natively on Ubuntu.  
> Windows users must use [WSL](https://learn.microsoft.com/en-us/windows/wsl/).

---

## Overview
SafeInCave is a 3D finite element simulator based on FEniCSx. It is designed to simulate the mechanical behavior of salt caverns under different operational conditions.

---

## Key Features

- **MPI-powered parallelism**: Scale simulations efficiently with mpi4py for distributed computing
- **Thermal effects**: Solve heat diffusion equation and include thermal strains and creep thermal responses
- **Pressure solution creep**: Easily add pressure solution creep to your simulations
- **Cavern abandonment**: Brine thermodynamic model and specific BC's for abandonment phase
- **Cyclic operations**: Impose fast cyclic pressure loads to the cavern walls
- **Constitutive model**: Include transient creep, reverse transient creep, and steady-state creep
- **Robust linearization**: Provides robustness and flexibility to include new constitutive models
- **Time discretization**: Choose between Explicit, Crank-Nicolson, and Fully-Implicit schemes
- **XDMF output**: Efficient output format in terms of size and postprocessing

---

## Getting started
The user should start by reading the documentation at folder docs/manual. There you will find the installation steps, simple to complex tutorials, and detailed instructions on how to setup your own simulation cases. In addition, check out our video lectures on the SafeInCave simulator:

1) Tensorial operations (theory): https://youtu.be/w5KX3F_rdzU?si=QQLVBq1NcrvOiS32
2) Tensorial operations (exercises): https://youtu.be/JiN6jwp0RPk?si=K1Qhe3lAxJD4LI5w
3) SafeInCave installation: https://youtu.be/J0KQ-nBMwwU?si=VXUVk60yrxavfT1E
4) Stay tuned for upcoming video lectures.

---

## Current members 
- [Hermínio Tasinafo Honório] (H.TasinafoHonorio@tudelft.nl),  Maintainer, 2023-present
- [Hadi Hajibeygi] (h.hajibeygi@tudelft.nl), Principal Investigator

---

## License
This project is licensed under the **GNU General Public License v3.0** (GPL-3.0).  
See the [LICENSE](LICENSE) file for full terms, or review the [official GPLv3 text](https://www.gnu.org/licenses/gpl-3.0.en.html).

---

## Papers and publications
[1] Honório, H.T, Houben, M., Bisdom, K., van der Linden, A., de Borst, K., Sluys, L.J., Hajibeygi, H. A multi-step calibration strategy for reliable parameter determination of salt rock mechanics constitutive models. Int J Rock Mech Min, 2024 (https://doi.org/10.1016/j.ijrmms.2024.105922)

[2] Honório, H.T, Hajibeygi, H. Three-dimensional multi-physics simulation and sensitivity analysis of cyclic hydrogen storage in salt caverns. Int J Hydrogen Energ, 2024 (https://doi.org/10.1016/j.ijhydene.2024.11.081)

[3] Kumar, K.R., Makhmutov, A., Spiers, C.J., Hajibeygi, H. Geomechanical simulation of energy storage in salt formations. Scientific Reports, 2022 (https://doi.org/10.1038/s41598-021-99161-8)

[4] Kumar, K.R., Hajibeygi, H. Influence of pressure solution and evaporate heterogeneity on the geo-mechanical behavior of salt caverns. The Mechanical Behavior of Salt X, 2022 (https://doi.org/10.1201/9781003295808)

## Acknowledgements
We would like to thank:
- Shell Global Solutions International B.V for sponsoring the project SafeInCave, within which this simulator was developed.
- Energi Simulation for currently supporting this project.

