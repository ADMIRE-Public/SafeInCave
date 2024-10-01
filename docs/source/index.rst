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


.. note::

   This project is under active development.

This simulator was developed as a deliverable of a project also called **SafeInCave** and financially supported by Shell.


Getting started
---------------

This chapter presents the installation steps for SafeInCave simulator and its dependencies. It also shows how to run a simple simulation of a triaxial test performed on a salt rock cubic sample.

.. toctree::

   getting_started



The Input File
---------------

The SafeInCave simulator runs entirely based on a single input file in JSON format. This chapter describes the structure of the input file and how it can be created automatically using class **InputFileAssistant**.

.. toctree::

   input_file



Tutorials
---------------

This chapter presents two tutorials that illustrate different capabilities of the SafeInCave simulator. The first tutorial handles a heterogeneous medium and shows how to assign different material properties for each region of the domain. The second tutorial addresses how to simulate a salt cavern with realistic boundary conditions.

.. toctree::

   tutorials


.. _fundamental-theory:

Fundamental Theory
------------------

This chapter is intended to provide the basic concepts of computational solid mechanics
necessary to understand the SafeInCave implementation.

.. toctree::

   introduction







API DOCUMENTATION
-----------------

API reference

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   modules