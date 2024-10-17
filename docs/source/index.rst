.. SafeInCave documentation master file, created by
   sphinx-quickstart on Sat Jul  6 11:37:00 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.

.

.. image:: _static/logo_2.png
   :alt: Project Logo
   :align: center
   :class: logo
   :width: 50%

Salt cavern simulations made easy!
==================================

The **SafeInCave** simulator is developed to study the mechanical behavior/stability of salt caverns for gas storage.
This is a three-dimensional simulator developed in Python, in which the finite element implementation is based on
`FEniCS 2019.1 <https://fenicsproject.org/download/archive/>`_.

**Features**:

- Fully open-source and documented.
- Comprehensive constitutive model able to capture transient creep, steady-state (dislocation) creep, and reverse transient creep.
- Robust numerical formulation for non-linear mechanics achieved by computing the consistent tangent matrix.
- Different choices of time integration schemes: explicit, Crank-Nicolson and fully-implicit.
- Interaction with the simulator happens through a single input file.
- No lines of code are necessary, except for building the input file (if necessary).
- Time-dependent and non-uniform boundary conditions can be assigned for the overburden, sideburden and gas pressure.
- Tests are added for the main classes and functions, thus enhancing robustness and facilitating external contributions.

**Current members**:

- [Hermínio Tasinafo Honório] (H.TasinafoHonorio@tudelft.nl),  Maintainer, 2023-present
- [Hadi Hajibeygi] (h.hajibeygi@tudelft.nl), Principal Investigator

**Acknowledgements**:

This simulator was developed as a deliverable of a project also called **SafeInCave** and financially supported by Shell.


Getting started
~~~~~~~~~~~~~~~~

This chapter presents the installation steps for SafeInCave simulator and its dependencies. It also shows how to run a simple simulation of a triaxial test performed on a salt rock cubic sample. See :ref:`tutorials` for more examples with detailed descriptions.

.. toctree::

   getting_started



The Input File
~~~~~~~~~~~~~~~~

The SafeInCave simulator runs entirely based on a single input file in JSON format. This chapter describes the structure of the input file and how it can be created automatically using class **InputFileAssistant**.

.. toctree::

   input_file


.. _tutorials:

Tutorials
~~~~~~~~~~~~~~~~

This chapter presents two tutorials that illustrate different capabilities of the SafeInCave simulator. The first tutorial handles a heterogeneous medium and shows how to assign different material properties for each region of the domain. The second tutorial addresses how to simulate a salt cavern with realistic boundary conditions.

.. toctree::

   tutorials


.. _fundamental-theory:

Fundamental Theory
~~~~~~~~~~~~~~~~~~

This chapter is intended to provide the basic concepts of computational solid mechanics
necessary to understand the SafeInCave implementation.

.. toctree::

   introduction







API DOCUMENTATION
~~~~~~~~~~~~~~~~~

API reference

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   modules