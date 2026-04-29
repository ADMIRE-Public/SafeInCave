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

## Install Conda + FEniCSx + SafeInCave

Use conda to avoid conflicts with other packages on your system. Download Miniconda3-py310_25.5.1-0-Linux-x86_64.sh from https://repo.anaconda.com/miniconda/ and save it in your Ubuntu home/user_name directory. In this same directory, execute:

```bash
bash Miniconda3-py310_25.3.1-1-Linux-x86_64.sh
```

Restart the terminal and execute:

```bash
conda activate
```

You should notice the tag (base) in your command line. Let's create a new environment named **safe**:

```bash
conda create -n safe python=3.10
conda activate safe
```

The tag (safe) should now appear in the command line. SafeInCave currently uses FEniCSx v0.9.0, so make sure to install this version. To install FEniCS v0.9.0 using conda (as described [here](https://fenicsproject.org/download/)), executed:

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

