# Installation

This section walks you through the installation process of the SafeInCave simulator on your system. Since SafeInCave is based on FEniCSx, it can only be installed on Ubuntu and MacOS systems. If you use Windows, first install Windows Subsystem for Linux (WSL) and then install Ubuntu on it. Next, install [FEniCSx](https://fenicsproject.org/), and finally the [safeincave](https://pypi.org/project/safeincave/) package.


## Install WSL + Ubuntu

**_NOTE:_** Skip this step if you are **not** using Windows.

To install WSL, open PowerShell in administrator mode and execute the following commands:

```bash
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
wsl --update
```

To proceed with Ubuntu installation, open PowerShell in administrator mode and execute:

```bash
wsl --install -d Ubuntu-22.04
```

Choose a username and a password. You will notice that PowerShell suddenly becomes a Ubuntu terminal.

## Install Conda

First, download Conda:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Then, install it by running the installer:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the prompt steps: 

- Press *Enter* to scroll
- Type *yes* to accept license
- Accept default install location (~/miniconda3 recommended)
- Say *yes* to `conda init`. 

Finally, close and reopen your WSL terminal.


## Install SafeInCave

Clone the SafeInCave repository to your local machine.

```bash
git clone https://github.com/ADMIRE-Public/SafeInCave.git
```

Go to *SafeInCave* root folder.

```bash
cd SafeInCave
```

It is recommended to use the latest stable release, currently 3.0.2.

```bash
git checkout v3.0.2
```

Install SafeInCave and all of its dependencies.

```bash
conda env create -f envs/user.yml
```

Activate the new conda environment.

```bash
conda activate safeincave-user
```

## Test installation

To test if installation was successful, you can run one of the examples in folder *examples*. For instance,

```bash
cd examples/mechanics/1_triaxial
python3 main.py
python3 plot_results.py
```



## Separate installation

Alternatively, the user may choose to manually install SafeInCave. The procedure is still simple, and it consists of installing FEniCSx 0.9.0 before installing SafeInCave with pip. 

First, install Conda (see Conda installation guidelines above).

Create a new environment named **safe**.

```bash
conda create -n safe python=3.10
conda activate safe
```

The tag (safe) should now appear in the command line. SafeInCave currently uses [FEniCSx v0.9.0](https://fenicsproject.org/blog/v0.9.0/), so make sure to install this version. To install FEniCS v0.9.0 using conda, executed:

```bash
conda install -c conda-forge fenics-dolfinx=0.9.0 mpich pyvista
```

Install SafeInCave by executing:

```bash
sudo apt update
sudo apt install python3-pip
pip install --upgrade pip
pip3 install safeincave
```

