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

It is advised to install Conda using Ubuntu terminal by following the instructions presented `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_. Once the Conda installation is finished, create the following environment on Ubuntu terminal:

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


Running your first simulation
------------------------------

The fastest way to run a simulation is to execute the *main.py* file of one of the examples in the *examples* folder. The example in folder *safeincave/examples/triaxial* simulates a triaxial test performed on a salt rock sample of cubic shape. The salt sample is subjected to a constant confining pressure and varying axial load, which can be visualized by executing the *plot_bcs.py* file, that is

.. code-block:: console

    (myenv) user@institution:~/safeincave$ cd examples/triaxial
    (myenv) user@institution:~/safeincave/examples/triaxial$ python plot_bcs.py

which produces the image shown in :numref:`Fig. %s <triaxial-load>`. 

.. _triaxial-load:

.. figure:: _static/triaxial_loading_schedule.png
   :align: center
   :width: 70%

   Loading schedule for the triaxial test example.

To run this example, simply do the following:

.. code-block:: console

    (myenv) user@institution:~/safeincave/examples/triaxial$ python main.py

Once the simulation is finished, the results can be found in folder *triaxial/output/case_0*. 

The results can be visualized on `Paraview <https://www.paraview.org/>`_. Alternatively, use `matplotlib <https://matplotlib.org/stable/>`_ to visualize results by doing the following:

.. code-block:: console

    (myenv) user@institution:~/safeincave/examples/triaxial$ python plot_results.py

This will generate :numref:`Fig. %s <triaxial-test>`, which shows the vertical (:math:`\varepsilon_v`) and horizontal (:math:`\varepsilon_h`) deformation over time.

.. _triaxial-test:

.. figure:: _static/triaxial_results.png
   :align: center
   :width: 80%

   Results obtained from the triaxial test simulation.

