import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from Grid import GridHandlerGMSH
from ResultsHandler import read_tensor_from_cells, read_tensor_from_cells_old
from Utils import *
from dolfin import *
import dolfin as do
import pandas as pd
import numpy as np

def save_dfs(case_folder):
	cells_coord, s_x, s_y, s_z, s_xy, s_xz, s_yz = read_tensor_from_cells_old(os.path.join(case_folder, "vtk", "stress"), "stress.pvd")

	if not os.path.exists(os.path.join(case_folder, "pandas")):
		os.makedirs(os.path.join(case_folder, "pandas"))

	cells_coord.to_pickle(os.path.join(case_folder, "pandas", "cells_coord.pkl"))
	s_x.to_pickle(os.path.join(case_folder, "pandas", "s_x.pkl"))
	s_y.to_pickle(os.path.join(case_folder, "pandas", "s_y.pkl"))
	s_z.to_pickle(os.path.join(case_folder, "pandas", "s_z.pkl"))
	s_xy.to_pickle(os.path.join(case_folder, "pandas", "s_xy.pkl"))
	s_xz.to_pickle(os.path.join(case_folder, "pandas", "s_xz.pkl"))
	s_yz.to_pickle(os.path.join(case_folder, "pandas", "s_yz.pkl"))


def exctract_material_props(input_model):
	n = np.array(input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["n"])[:,np.newaxis]
	gamma = np.array(input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["gamma"])[:,np.newaxis]
	beta_1 = np.array(input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["beta_1"])[:,np.newaxis]
	beta = np.array(input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["beta"])[:,np.newaxis]
	m = np.array(input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["m"])[:,np.newaxis]
	sigma_t = np.array(input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["sigma_t"])[:,np.newaxis]
	return n, gamma, beta_1, beta, m, sigma_t

def read_stresses(case_folder):
	s_xx = pd.read_pickle(os.path.join(case_folder, "s_x.pkl"))
	s_xx = pd.read_pickle(os.path.join(case_folder, "s_x.pkl"))
	s_yy = pd.read_pickle(os.path.join(case_folder, "s_y.pkl"))
	s_zz = pd.read_pickle(os.path.join(case_folder, "s_z.pkl"))
	s_xy = pd.read_pickle(os.path.join(case_folder, "s_xy.pkl"))
	s_xz = pd.read_pickle(os.path.join(case_folder, "s_xz.pkl"))
	s_yz = pd.read_pickle(os.path.join(case_folder, "s_yz.pkl"))
	return s_xx, s_yy, s_zz, s_xy, s_xz, s_yz

def compute_FOS(case_folder):
	try:
		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = read_stresses(os.path.join(case_folder, "pandas"))
	except:
		save_dfs(case_folder)
		s_xx, s_yy, s_zz, s_xy, s_xz, s_yz = read_stresses(os.path.join(case_folder, "pandas"))

	# Extract material properties
	input_model = read_json(os.path.join(case_folder, "input_file.json"))
	n, gamma, beta_1, beta, m, sigma_t = exctract_material_props(input_model)

	# Extract time steps
	times = s_xx.columns.values

	# Build numpy array stresses in MPa
	s_xx = -s_xx.values/MPa
	s_yy = -s_yy.values/MPa
	s_zz = -s_zz.values/MPa
	s_xy = -s_xy.values/MPa
	s_xz = -s_xz.values/MPa
	s_yz = -s_yz.values/MPa

	# Compute stress invariants
	I1 = s_xx + s_yy + s_zz
	I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
	I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
	J2 = (1/3)*I1**2 - I2
	J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
	Sr = -(J3*np.sqrt(27))/(2*J2**1.5)

	# Compute Factor Of Safety (FOS)
	F_dil = (1 - 2/n)*(gamma*I1**2)*(np.exp(beta_1*I1) + beta*Sr)**m
	FOS = np.sqrt(F_dil)/np.sqrt(J2)
	return times, FOS


def main():
	# Define case folder
	case_folder = os.path.join("output", "case_0", "operation")

	# Load mesh
	mesh_folder = os.path.join("..", "..", "grids", "cavern_irregular")
	g = GridHandlerGMSH("geom", mesh_folder)

	P0 = do.FunctionSpace(g.mesh, "DG", 0)
	FOS = do.Function(P0)
	FOS.rename("Factor of Safety", "-")

	FOS_vtk = do.File(os.path.join(case_folder, "vtk", "FOS", "FOS.pvd"))

	times, FOS_data = compute_FOS(case_folder)

	for i, time in enumerate(times):
		FOS.vector()[:] = FOS_data[:,i]
		FOS_vtk << (FOS, time)


if __name__ == '__main__':
	main()