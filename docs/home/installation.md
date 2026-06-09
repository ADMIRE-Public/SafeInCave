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

Follow the prompt steps: Press *Enter* to scroll, type *yes* to accept license, accept default install location (~/miniconda3 recommended), and say *yes* to `conda init`. 

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
