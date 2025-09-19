<img src="assets/logo_2.png" width="500"/>

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
- **Graphical user interface**: Build your simulation without writing lines of code
- **Constitutive model**: Include transient, reverse transient, dislocation, and pressure solution creep
- **Robust linearization**: Provides robustness and flexibility to include new constitutive models
- **Time discretization**: Choose between Explicit, Crank-Nicolson, and Fully-Implicit schemes
- **XDMF output**: Efficient output format in terms of size and postprocessing

---

## Installation
SafeInCave installation depends on [FEniCSx](https://fenicsproject.org/) installation. For Windows users, the installaion pipeline consists of:

1) Install [WSL](https://learn.microsoft.com/en-us/windows/wsl/)

2) Install Ubuntu

3) Install [FEniCSx](https://fenicsproject.org/download/)

4) Install SafeInCave

See SafeInCave [documentation](https://safeincave-docs.streamlit.app/installation) for a detailed explanation on the installation process.

---

## Getting started
After installation, the easiest way to set up SafeInCave simulations is by using the SafeInCave App, as shown in the image below.

<img src="assets/gui_safeincave.jpeg" width="800"/>

Alternatively, users can build their own simulators using the *safeincave* package. Our [documentation](https://safeincave-docs.streamlit.app/) webpage shows detailed examples of how to set up purely [mechanical]() simulations, [heat diffusion]() simulations, and [thermomechanical]() simulations. These examples show how to build constitutive models, apply different types of boundary conditions, assign material properties, etc.

---

## Extra material
Video lectures and video tutorials can be found in the [ADMIRE](https://www.youtube.com/@ADMIRE1/featured) YouTube channel. The following videos are currently available:

1) [Tensorial operations (theory)](https://youtu.be/w5KX3F_rdzU?si=QQLVBq1NcrvOiS32)

2) [Tensorial operations (exercises)](https://www.youtube.com/watch?v=JiN6jwp0RPk&t=0s)

3) [Constitutive modeling](https://www.youtube.com/watch?v=fCeJIbjIL10)

4) Stay tuned to [ADMIRE](https://www.youtube.com/@ADMIRE1/featured) YouTube channel for upcoming video lectures.

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

---

## Acknowledgements
We would like to thank:
- Shell Global Solutions International B.V for sponsoring the project SafeInCave, within which this simulator was developed.
- Energi Simulation for currently supporting this project.

