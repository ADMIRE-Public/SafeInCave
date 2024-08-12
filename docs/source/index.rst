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

The **SafeInCave** simulator is intended to study the mechanical behavior/stability of salt caverns for gas storage.
This is a three-dimensional simulator developed in Python, in which the finite element implementation is based on
`FEniCS 2019.1 <https://fenicsproject.org/download/archive/>`_.


Check out the :ref:`fundamental-theory` section for further information.

.. note::

   This project is under active development.



Getting started
---------------

.. toctree::

   getting_started



Tutorials
---------------

This chapter presents some tutorials.

.. admonition:: Notation:

    - Position vector: :math:`\mathbf{r} = \begin{bmatrix} x & y & z \end{bmatrix}^T`.
    - Displacement vector: :math:`\mathbf{u} = \begin{bmatrix} u & v & w \end{bmatrix}^T`.
    - Normal vector: :math:`\mathbf{n} = \begin{bmatrix} n_x & n_y & n_z \end{bmatrix}^T`.
    - Stress tensor: :math:`\pmb{\sigma} = \begin{bmatrix} \sigma_{xx} & \sigma_{xy} & \sigma_{xz} \\ \sigma_{xy} & \sigma_{yy} & \sigma_{yz} \\ \sigma_{xz} & \sigma_{yz} & \sigma_{zz} \end{bmatrix}^T`.

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