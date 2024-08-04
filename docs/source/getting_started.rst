

Installation
------------

The SafeInCave simulator requires the following dependencies to be installed:

- `PyTorch <https://pytorch.org/>`_
- `pandas <https://pandas.pydata.org/>`_
- `FEniCS 2019.1 <https://fenicsproject.org/download/archive/>`_

The SafeInCave simulator has been developed and tested only on Windows platform, but it should also work with other operational systems. In this section, only Windows installation is covered.

The SafeInCave simulator is based on `FEniCS 2019.1 <https://fenicsproject.org/download/archive/>`_, 

Because the SafeInCave simulator is based on `FEniCS 2019.1 <https://fenicsproject.org/download/archive/>`_, which can be installed on Ubuntu, the installation on Windows requires the Windows Subsystem for Linux (WSL). To install WSL, open the Power Shell in **administrator mode** and run the following commands:

.. code-block:: powershell

    Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
    wsl --update

To install Ubuntu, run the following command on Power Shell:

.. code-block:: powershell

    wsl --install -d ubuntu

At this point, the Ubuntu terminal should be available at *Start button* -> *Ubuntu*. 

It is advised to install Conda using Ubuntu terminal by following the instructions presented `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_. When the Conda installation is finished, create the following environment on Ubuntu terminal:

.. code-block:: console

    (base) $ conda create --name myenv python=3.8
    (base) $ conda activate myenv
    (myenv) $ python --version
    Python 3.8

To install PyTorch, follow these `instructions <https://pytorch.org/>`_.

To install pandas, run the following commands on Ubuntu terminal:

.. code-block:: console

    (myenv) $ conda install pandas
    (myenv) $ conda install -c anaconda openpyxl

To install FEniCS 2019.1, run the following command:

.. code-block:: console

    (myenv) $ conda install -c conda-forge fenics=2019.1

Download `Gmsh version 4.10.1 <https://gmsh.info/bin/Windows/>`_




Input File
----------
The SafeInCave simulator runs entirely based on a single input file in JSON format. This file can be either created manually or with the help of class **InputFileAssistant**. The input file requires the following sections:

1. *grid*
2. *output*
3. *solver_settings*
4. *simulation_settings*
5. *body_force*
6. *time_settings*
7. *boundary_conditions*
8. *constitutive_model*

Therefore, the input file should have the basic structure shown below.

.. code-block:: json

    {
        "grid": {},
        "output": {},
        "solver_settings": {},
        "time_settings": {},
        "simulation_settings": {},
        "body_force": {},
        "boundary_conditions": {},
        "constitutive_model": {}
    }

A detailed explanation of each section is presented next.

Section *grid*
~~~~~~~~~~~~~~
The section grid informs the grid to be used in the simulation. This section requires two keys: (1) *path* and (2) *name*. The key path specifies the relative path to the directory where the grid stored. The key *name* indicates the name of the grid files, which is usually *geom* (e.g. *geom.xml*, *geom_facet_region.msh*, *geom_physical_region.xml*). The snippet :numref:`Listing %s <grid-section>` illustrates a typical example.

.. _grid-section:

.. code-block:: json
    :caption: Input file section: *grid*.

    {
        "grid": {
            "path": "../../grids/cube_0",
            "name": "geom"
        },
    }

Section *output*
~~~~~~~~~~~~~~~~
The section *output* only requires the key *path*, which is specifies the relative path to the directory where the output files will be saved. This is illustrated in :numref:`Listing %s <output-section>`, where the results are saved in the directory *output/case_name*. This directory is automatically created in case it does not exist.

.. _output-section:

.. code-block:: json
    :caption: Input file section: *output*

    {
        "output": {
            "path": "output/case_name"
         },
    }

Section *solver_settings*
~~~~~~~~~~~~~~~~~~~~~~~~~
This section specifies which solver is used to solve the linear systems. The required keys for *solver_settings* are *type* and *method*. The *type* key can be either *LU* or *KrylovSolver*, for direct LU decomposition of a Krylov-based solver, respectively. If *LU* is chosen, the *method* key can be either *default*, *umfpack*, *mumps*, *pastix*, *superlu*, *superlu_dist*, or *petsc*, dependint on how PETSc has been installed. For example,

.. _solver-settings-lu:

.. code-block:: json
    :caption: Input file section: *solver_settings* (LUSolver)

    {
        "solver_settings": {
            "type": "LU",
            "method": "petsc"
         },
    }

A Krylov-based solver can be chosen by specifying the keyword *KrylovSolver* to the *type* key. The specific Krylov solver is defined under the key *method*, and the main options are: *cg*, *bicg*, *bigcstab*, and *gmres*. In addition to *type* and *method*, the *KrylovSolver* requires keys *preconditioner* and *relative_tolerance*. The main options for key *preconditioner* are: *icc*, *ilu*, *petsc_amg*, *sor*, and *hypre*. For example,

.. _solver-settings-krylov:

.. code-block:: json
    :caption: Input file section: *solver_settings* (KrylovSolver)

    {
        "solver_settings": {
            "type": "KrylovSolver",
            "method": "cg",
            "preconditioner": "petsc_amg",
            "relative_tolerance": 1e-12
         },
    }

Section *simulation_settings*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This section specifies whether or not to compute the equilibrium condition before the actual simulation begins. It requires the *equilibrium* and *operation* keys, which specifies the settings for the equilibrium and operation simulation stages, respectively. In the equilibrium condition, the stresses specified at the initial time :math:`t=0` are applied to the geometry and a simulation is run considering only the **elastic** and **viscoelastic** (if present) part of the constitutive model. This equilibrium simulation is run until the it reaches steady-state condition. The keyword *true* or *false* specify whether the equilibrium condition is computed or not. The key *dt_max* specifies the time step size adopted to reach steady-state condition, which is defined by the *time_tol* key.

The *operation* key requires the key *active*, which can be *true* or *false*. The *dt_max* key defines the time step size of the simulation during the operation stage. Finally the *n_skip* key specifies how many time steps to skip before saving the results. This is useful in simulations where a very small time step size is required, thus avoiding excessively large results files. An example is shown in :numref:`Listing %s <simulation-settings>`.

.. _simulation-settings:

.. code-block:: json
    :caption: Input file section: *simulation_settings*

    {
        "simulation_settings": {
           "equilibrium": {
               "active": true,
               "dt_max": 1800.0,
               "time_tol": 0.0001
           },
           "operation": {
               "active": true,
               "dt_max": 1800.0,
               "n_skip": 1
           }
        },
    }

Section *body_forces*
~~~~~~~~~~~~~~~~~~~~~
This section defines the body forces associated to the rock mass. The gravity acceleration is specified under the key *gravity*; the rock density is defined under key *ensity*; and the direction along which the gravity acceleration is aligned is specified under the key *direction* (0 for *x*, 1 for *y* and 2 for *z*). For example, :numref:`Listing %s <body-force>`.

.. _body-force:

.. code-block:: json
    :caption: Input file section: *body_force*

    {
        "body_force": {
            "gravity": -9.81,
            "density": 2000,
            "direction": 2
        },
    }

.. _time-settings:

Section *time_settings*
~~~~~~~~~~~~~~~~~~~~~~~
In the *time_settings* section, the time integration method is defined by chosing the :math:`\theta` value under the key *theta* (0 for fully-implicit, 0.5 for Crank-Nicolson, and 1 for explicit). Next, the key *timeList* specifies the time schedule that defines the loading conditions (see :ref:`section-boundary-conditions`). For example,

.. _time-settings-section:

.. code-block:: json
    :caption: Input file section: *time_settings*

    {
        "time_settings": {
            "theta": 0.0,
            "time_list": [0, 10, 20]
        },
    }

.. _section-boundary-conditions:

Section *boundary_conditions*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This section allows for specifying the boundary conditions of the problem. For salt cavern simulations, it is often the case that the pressure inside the cavern varies with time. Additionally, for very tall caverns, there is a significant pressure difference between the top and the bottom of the cavern due to the gas specific weight. The sideburden, althought fixed in time, also varies significantly from top to bottom of the geometry. The section *boundary_conditions* was designed to allow for an easy spefication of such boundary conditions. To exemplify this process, consider the examples illustrated in :numref:`Fig. %s <bc-block-full>`, which shows a 2D view of a block with boundaries names *BOTTOM*, *TOP*, *WEST* and *EAST*. :numref:`Fig. %s <bc-block-full>`-a shows in details the boundary conditions applied at the initial time step :math:`t_0`. As it can be verified, the *BOTTOM* and *WEST* boundaries are prevented from normal displacement (Dirichlet boundary condition), whereas the *TOP* boundary is subjected to a constant (in space) compressive load, and a *z*-dependent load is applied to boundary *EAST*. Moreover, :numref:`Fig. %s <bc-block-full>`-b shows that the applied loads actually vary with time.

.. _bc-block-full:

.. figure:: _static/bc_block_full.png
   :alt: block
   :align: center
   :width: 100%

   Boundary conditions applied to block.

The keys inside the *boundary_settings* section must be the boundary names. Inside each boundary name, there is a *type* key that can be either *dirichlet* or *neumann*. If *type* is *dirichlet, then the imposed displacement component must be specified under the key *component* (0 for *x*, 1 for *y* and 2 for *z*). Next, the key *values* receives a list of prescribed values for each time level according to the *time_list*, defined in section *time_settings* (**both lists must be the same size**). If *type* is *neumann*, then the keys *direction*, *density*, *reference_position* and *values* are required. The *direction* key defines the direction along which the boundary condition varies spacially; the *density* key specifies how much the load changes in that direction; the *reference_position* key defines the position :math:`H` where the specified values :math:`p_0` are located (see :numref:`Fig. %s <bc-block-full>`-a); and the *values* key receives a list of prescribed loads corresponding to each time of *time_settings*.

The boundary conditions illustrated in :numref:`Fig. %s <bc-block-full>` are written in the JSON file as shown below (:numref:`Listing %s <boundary-conditions>`). The *BOTTOM* and *WEST* boundaries are of type *dirichlet* with value 0 in the time interval between 0 and 20 s (see :ref:`time-settings`). The displacement component normal to boundary *BOTTOM* is in the *z* direction, that is why the key *component* receives the value 2. On the other hand, the normal displacement on boundary *WEST* is aligned to the *x* direction, thus the value 0 to the key *component*. The boundary *EAST* is subjected to a boundary condition of *type* *neumann*, and the spatial variation takes place in the *z* direction (*direction: 2*). The amount of variation :math:`\rho` is specified as *density: 50* and the reference position :math:`H` is *reference_position: 1.0*, according to :numref:`Fig. %s <bc-block-full>`-a. According to :numref:`Fig. %s <bc-block-full>`-a, the load imposed on the *TOP* boundary is uniform, so the *density* key should be zero. As a consequence, the value specified in the *direction* and *reference_position* keys and do not matter at all.

.. note::
    
    The value of gravity :math:`g` shown in :numref:`Fig. %s <bc-block-full>`-a is specified in :ref:`body-force`.

.. _boundary-conditions:

.. code-block:: json
    :caption: Input file section: *boundary_conditions*

    {
        "boundary_conditions": {
            "BOTTOM": {
                "type": "dirichlet",
                "component": 2,
                "values": [0.0, 0.0, 0.0]
            },
            "WEST": {
                "type": "dirichlet",
                "component": 0,
                "values": [0.0, 0.0, 0.0]
            },
            "EAST": {
                "type": "neumann",
                "density": 50.0,
                "direction": 2,
                "reference_position": 1.0,
                "values": [5.0, 7.0, 10.0]
            },
            "TOP": {
                "type": "neumann",
                "density": 0.0,
                "direction": 0,
                "reference_position": 0.0,
                "values": [5.0, 8.0, 5.0]
            }
        }
    }


Section *constitutive_model*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SafeInCave simulator allows for very flexible choices of the constitutive model. As an example, we consider the constitutive model illustrated in :numref:`Fig. %s <constitutive-model-0>`, which is composed of a linear spring element, two Kelvin-Voigt elements, one viscoplastic element, and one dislocation creep element. Each one of these elements comprise its own set of material parameters, as indicated in the figure. Refer to :ref:`constitutive-models-section` for a detailed explanation of each element and the corresponding material properties.

.. _constitutive-model-0:

.. figure:: _static/constitutive_model_0.png
   :alt: block
   :align: center
   :width: 75%

   Elements composing the constitutive model.

To be general, let us consider a simple mesh divided in two sub-domains with different material properties. This is illustrated in :numref:`Fig. %s <mesh-regions>`, where elements 0, 1, 4 and 5 belong to :math:`\Omega_A`, while elements 2, 3, 6 and 7 belong to :math:`\Omega_B`.

.. note::

    A 2D grid is considered here only for simplicity. However, the SafeInCave simulator only handles 3D grids composed of tetrahedral elements.


.. _mesh-regions:

.. figure:: _static/mesh_regions.png
   :alt: block
   :align: center
   :width: 35%

   Computational mesh divided in two sub-domains: :math:`\Omega_A` and :math:`\Omega_B`.

The material properties assigned to each sub-domain is presented in :numref:`Table %s <table-mat-props>`. In this example, the values assigned to each material property are merely illustrative and **do not** correspond real physical values.

.. _table-mat-props:

.. list-table:: Material properties for domains :math:`\Omega_A` and :math:`\Omega_B`.
   :widths: 25 25 25
   :header-rows: 1

   * - Property name
     - Domain :math:`\Omega_A`
     - Domain :math:`\Omega_B`
   * - :math:`E_0`
     - 100
     - 250
   * - :math:`\nu_0`
     - 0.3
     - 0.2
   * - :math:`E_1`
     - 90
     - 75
   * - :math:`\nu_1`
     - 0.15
     - 0.42
   * - :math:`\eta_1`
     - 7.0
     - 8.2
   * - :math:`E_2`
     - 120
     - 165
   * - :math:`\nu_2`
     - 0.24
     - 0.38
   * - :math:`\eta_2`
     - 17.0
     - 6.3
   * - :math:`\mu_1`
     - 5.3
     - 2.1
   * - :math:`N_1`
     - 3.1
     - 2.9
   * - :math:`n`
     - 3
     - 3
   * - :math:`a_1`
     - 1.9
     - 2.3
   * - :math:`\eta`
     - 0.82
     - 0.97
   * - :math:`\beta`
     - 0.99
     - 0.76
   * - :math:`\beta_1`
     - 0.38
     - 0.75
   * - :math:`m`
     - -0.5
     - -0.5
   * - :math:`\gamma`
     - 0.087
     - 0.095
   * - :math:`\alpha_0`
     - 0.40
     - 0.27
   * - :math:`k_v`
     - 0.0
     - 0.6
   * - :math:`\sigma_t`
     - 5.0
     - 4.2
   * - :math:`A_1`
     - 1.9
     - 2.3
   * - :math:`n_1`
     - 3.1
     - 4.2
   * - :math:`T`
     - 298
     - 298
   * - :math:`Q`
     - 51600
     - 51600
   * - :math:`R`
     - 8.32
     - 8.32



The *constitutive_model* section requires three mandatory keys: *Elastic*, *Viscoelastic* and *Inelastic*. A spring can be added to the *Elastic* key as shown in :numref:`Listing %s <constitutive-model>`. The name *Spring0* is an arbitrary name given to the spring; the *type* key must be *Spring*; the key *active* can be *true* or *false* depending on whether the user wants to include it or not to the constitutive model; finally, the key *parameters* contains the lists of the material parameters associated to the spring (i.e. Young's modulus, :math:`E`, and Poisson's ratio, :math:`\nu`). The size of these lists must be the same as the number of grid elements (in this case, 8 elements, as shown in :numref:`Fig. %s <mesh-regions>`). Therefore, the values in these lists represent the material properties of each element of the grid.

.. important::

    A constitutive model **must** include at least one spring. In other words, at least one spring must be **active**.

A Kelvin-Voigt element is a parallel arrangement between a linear spring and a linear dashpot. This type of element is added under the key *Viscoelastic*. In the example shown in :numref:`Listing %s <constitutive-model>`, two Kelvin-Voigt elements are added, namely, *KelvinVoigt1* and *KelvinVoigt2*. The key *type* must be *KelvinVoigt*. The material parameters associated to the Kelvin-Voigt element are the Poisson's ratio (:math:`\nu`) and Young's modulus (:math:`E`) of the spring, and the viscosity (:math:`\eta`) of the dashpot.

.. note::

    A Kelvin-Voigt element with a nonlinear dashpot, if implemented, should be also added under the *Viscoelastic* key.

The viscoplastic and dislocation creep elements in :numref:`Fig. %s <constitutive-model-0>` must be included under the *Inelastic* key. In this example, the arbitrary names given to the viscoplastic and dislocation creep elements are *ViscPlastDesai* and *DisCreep*, respectively. The viscoplastic element must be of type *ViscoplasticDesai* and dislocation creep element must be of type *DislocationCreep*.



.. _constitutive-model:

.. code-block:: json
    :caption: Input file section: *constitutive_model*

    {
        "constitutive_model": {
            "Elastic": {
                "Spring0": {
                    "type": "Spring",
                    "active": true,
                    "parameters": {
                        "E": [100, 100, 250, 250, 100, 100, 250, 250],
                        "nu": [0.3, 0.3, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2]
                    }
                }
            },
            "Viscoelastic": {
                "KelvinVoigt1": {
                    "type": "KelvinVoigt",
                    "active": true,
                    "parameters": {
                        "E": [90.0, 90.0, 75.0, 75.0, 90.0, 90.0, 75.0, 75.0],
                        "nu": [0.15, 0.15, 0.42, 0.42, 0.15, 0.15, 0.42, 0.42],
                        "eta": [7.0, 7.0, 8.2, 8.2, 7.0, 7.0, 8.2, 8.2]
                    }
                },
                "KelvinVoigt2": {
                    "type": "KelvinVoigt",
                    "active": true,
                    "parameters": {
                        "E": [120.0, 120.0, 165.0, 165.0, 120.0, 120.0, 165.0, 165.0],
                        "nu": [0.24, 0.24, 0.38, 0.38, 0.24, 0.24, 0.38, 0.38],
                        "eta": [17.0, 17.0, 6.3, 6.3, 17.0, 17.0, 6.3, 6.3]
                    }
                }
            },
            "Inelastic": {
                "ViscPlastDesai": {
                    "type": "ViscoplasticDesai",
                    "active": true,
                    "parameters": {
                        "mu_1": [5.3, 5.3, 2.1, 2.1, 5.3, 5.3, 2.1, 2.1],
                        "N_1": [3.1, 3.1, 2.9, 2.9, 3.1, 3.1, 2.9, 2.9],
                        "n": [3, 3, 3, 3, 3, 3, 3, 3],
                        "a_1": [1.9, 1.9, 2.3, 2.3, 1.9, 1.9, 2.3, 2.3],
                        "eta": [0.82, 0.82, 0.97, 0.97, 0.82, 0.82, 0.97, 0.97],
                        "beta_1": [0.99, 0.99, 0.76, 0.76, 0.99, 0.99, 0.76, 0.76],
                        "beta": [0.38, 0.38, 0.75, 0.75, 0.38, 0.38, 0.75, 0.75],
                        "m": [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5],
                        "gamma": [0.087, 0.087, 0.095, 0.095, 0.087, 0.087, 0.095, 0.095],
                        "alpha_0": [0.40, 0.40, 0.27, 0.27, 0.40, 0.40, 0.27, 0.27],
                        "k_v": [0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 0.6, 0.6],
                        "sigma_t": [5.0, 5.0, 4.2, 4.2, 5.0, 5.0, 4.2, 4.2]
                    }
                },
                "DisCreep": {
                    "type": "DislocationCreep",
                    "active": true,
                    "parameters": {
                        "A": [1.9, 1.9, 2.3, 2.3, 1.9, 1.9, 2.3, 2.3],
                        "n": [3.1, 3.1, 4.2, 4.2, 3.1, 3.1, 4.2, 4.2],
                        "T": [298, 298, 298,298, 298, 298, 298,298],
                        "Q": [51600, 51600, 51600, 51600, 51600, 51600, 51600, 51600],
                        "R": [832, 832, 832, 832, 832, 832, 832, 832]
                    }
                }
            }
        }
    }

The elements available for composing the constitutive model are summarized in :numref:`Table %s <list-elements>`, where the correspoding material parameters are also shown. The parameters, as discussed above, must be informed as a list of values associated to each grid element. Currently, the linear elastic spring, the viscoelastic Kelvin-Voigt element, the viscoplastic model of Desai (1987), and the dislocation creep element are implemented in the SafeInCave simulator.

.. _list-elements:

.. list-table:: Available elements for the constitutive model.
   :widths: 5 8 25
   :header-rows: 1

   * - Category
     - Type
     - Material parameters
   * - Elastic
     - Spring
     - E, nu
   * - Viscoelastic
     - KelvinVoigt
     - E, nu, eta
   * - Inelastic
     - DislocationCreep
     - A, n, T, Q, R
   * - Inelastic
     - ViscoplasticDesai
     - mu_1, N_1, n, a_1, eta, beta_1, beta, m, gamma, alpha_0, k_v, sigma_t

Full input file
~~~~~~~~~~~~~~~

To conclude this section, the complete input file should look like as in :numref:`Listing %s <full-input-file>`

.. _full-input-file:

.. code-block:: json
    :caption: Complete input file

    {
        "grid": {
            "path": "../../grids/cube_0",
            "name": "geom"
        },
        "output": {
            "path": "output/case_name"
         },
        "solver_settings": {
            "type": "KrylovSolver",
            "method": "cg",
            "preconditioner": "petsc_amg",
            "relative_tolerance": 1e-12
         },
        "simulation_settings": {
            "equilibrium": {
               "active": true,
               "dt_max": 1800.0,
               "time_tol": 0.0001
            },
            "operation": {
               "active": true,
               "dt_max": 1800.0,
               "n_skip": 1
            }
        },
        "body_force": {
            "gravity": -9.81,
            "density": 2000,
            "direction": 2
        },
        "time_settings": {
            "theta": 0.0,
            "time_list": [0, 10, 20]
        },
        "boundary_conditions": {
            "BOTTOM": {
                "type": "dirichlet",
                "component": 2,
                "values": [0.0, 0.0, 0.0]
            },
            "WEST": {
                "type": "dirichlet",
                "component": 0,
                "values": [0.0, 0.0, 0.0]
            },
            "EAST": {
                "type": "neumann",
                "density": 50.0,
                "direction": 2,
                "reference_position": 1.0,
                "values": [5.0, 7.0, 10.0]
            },
            "TOP": {
                "type": "neumann",
                "density": 0.0,
                "direction": 0,
                "reference_position": 0.0,
                "values": [5.0, 8.0, 5.0]
            }
        },
        "constitutive_model": {
            "Elastic": {
                "Spring0": {
                    "type": "Spring",
                    "active": true,
                    "parameters": {
                        "E": [100, 100, 250, 250, 100, 100, 250, 250],
                        "nu": [0.3, 0.3, 0.2, 0.2, 0.3, 0.3, 0.2, 0.2]
                    }
                }
            },
            "Viscoelastic": {
                "KelvinVoigt1": {
                    "type": "KelvinVoigt",
                    "active": true,
                    "parameters": {
                        "E": [90.0, 90.0, 75.0, 75.0, 90.0, 90.0, 75.0, 75.0],
                        "nu": [0.15, 0.15, 0.42, 0.42, 0.15, 0.15, 0.42, 0.42],
                        "eta": [7.0, 7.0, 8.2, 8.2, 7.0, 7.0, 8.2, 8.2]
                    }
                },
                "KelvinVoigt2": {
                    "type": "KelvinVoigt",
                    "active": true,
                    "parameters": {
                        "E": [120.0, 120.0, 165.0, 165.0, 120.0, 120.0, 165.0, 165.0],
                        "nu": [0.24, 0.24, 0.38, 0.38, 0.24, 0.24, 0.38, 0.38],
                        "eta": [17.0, 17.0, 6.3, 6.3, 17.0, 17.0, 6.3, 6.3]
                    }
                }
            },
            "Inelastic": {
                "ViscPlastDesai": {
                    "type": "ViscoplasticDesai",
                    "active": true,
                    "parameters": {
                        "mu_1": [5.3, 5.3, 2.1, 2.1, 5.3, 5.3, 2.1, 2.1],
                        "N_1": [3.1, 3.1, 2.9, 2.9, 3.1, 3.1, 2.9, 2.9],
                        "n": [3, 3, 3, 3, 3, 3, 3, 3],
                        "a_1": [1.9, 1.9, 2.3, 2.3, 1.9, 1.9, 2.3, 2.3],
                        "eta": [0.82, 0.82, 0.97, 0.97, 0.82, 0.82, 0.97, 0.97],
                        "beta_1": [0.99, 0.99, 0.76, 0.76, 0.99, 0.99, 0.76, 0.76],
                        "beta": [0.38, 0.38, 0.75, 0.75, 0.38, 0.38, 0.75, 0.75],
                        "m": [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5],
                        "gamma": [0.087, 0.087, 0.095, 0.095, 0.087, 0.087, 0.095, 0.095],
                        "alpha_0": [0.40, 0.40, 0.27, 0.27, 0.40, 0.40, 0.27, 0.27],
                        "k_v": [0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 0.6, 0.6],
                        "sigma_t": [5.0, 5.0, 4.2, 4.2, 5.0, 5.0, 4.2, 4.2]
                    }
                },
                "DisCreep": {
                    "type": "DislocationCreep",
                    "active": true,
                    "parameters": {
                        "A": [1.9, 1.9, 2.3, 2.3, 1.9, 1.9, 2.3, 2.3],
                        "n": [3.1, 3.1, 4.2, 4.2, 3.1, 3.1, 4.2, 4.2],
                        "T": [298, 298, 298,298, 298, 298, 298,298],
                        "Q": [51600, 51600, 51600, 51600, 51600, 51600, 51600, 51600],
                        "R": [832, 832, 832, 832, 832, 832, 832, 832]
                    }
                }
            }
        }
    }

Class **InputFileAssistant**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Explain the **InputFileAssistant** class.



Example 1
~~~~~~~~~