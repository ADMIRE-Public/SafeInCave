import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from Grid import GridHandlerGMSH
from ResultsHandler import read_tensor_from_cells, read_tensor_from_cells_old, read_scalar_from_cells
from Utils import *
from dolfin import *
import dolfin as do
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt



dilation_points = np.array([
				[-5.289256198347102, 0.06089309878213811],
				[-3.3057851239669382, 1.8876860622462814],
				[-0.11019283746556141, 3.897158322056839],
				[3.4159779614325068, 5.7848443843031205],
				[6.611570247933887, 7.4289580514208495],
				[10.24793388429752, 8.951285520974295],
				[13.553719008264462, 10.351826792963472],
				[17.190082644628095, 11.81326116373478],
				[20.495867768595044, 13.152909336941818],
				[23.80165289256198, 14.431664411366718],
				[27.21763085399449, 15.588633288227339],
				[30.523415977961424, 16.745602165087963],
				[34.0495867768595, 17.96346414073072],
				[37.465564738292, 19.120433017591342],
				[40.881542699724505, 20.15561569688769],
				[44.29752066115701, 21.190798376184034],
				[47.82369146005509, 22.34776725304466],
				[51.23966942148759, 23.322056833558864],
				[54.87603305785122, 24.41813261163735],
				[58.402203856749296, 25.39242219215156],
				[61.8181818181818, 26.305818673883632],
				[65.45454545454544, 27.341001353179976],
				[68.9807162534435, 28.254397834912048],
				[72.50688705234158, 29.289580514208392],
				[76.03305785123966, 30.202976995940464],
				[79.33884297520659, 30.994587280108256] ])
p_points = dilation_points[:,0]/3
q_points = dilation_points[:,1]*np.sqrt(3)


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
	n = input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["n"]
	gamma = input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["gamma"]
	beta_1 = input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["beta_1"]
	beta = input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["beta"]
	m = input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["m"]
	sigma_t = input_model["constitutive_model"]["inelastic"]["desai"]["parameters"]["sigma_t"]
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
	cells_coord, p = read_scalar_from_cells(os.path.join(case_folder, "vtk", "p"), "p.pvd")
	cells_coord, q = read_scalar_from_cells(os.path.join(case_folder, "vtk", "q"), "q.pvd")

	times = q.columns.values

	p = -p.values
	q = q.values

	FOS = np.zeros_like(q)
	for time_id in range(len(times)):
		for cell_id in range(q.shape[0]):
			q_dil = np.interp(p[cell_id, time_id], p_points, q_points)
			FOS[cell_id, time_id] = q_dil / q[cell_id, time_id]

	
	# fig, ax = plt.subplots(1, 1, figsize=(5, 3))
	# fig.subplots_adjust(top=0.970, bottom=0.155, left=0.09, right=0.980, hspace=0.35, wspace=0.225)

	# ax.plot(p_points, q_points, "k-")
	# ax.grid(True)

	# time_id = 10
	# ax.scatter(p[:,time_id], q[:,time_id], marker=".")

	return times, FOS


def main():
	# Define case folder
	case_folder = os.path.join("output", "case_0", "operation")

	# Load mesh
	mesh_folder = os.path.join("..", "..", "grids", "cavern_overburden")
	g = GridHandlerGMSH("geom", mesh_folder)

	P0 = do.FunctionSpace(g.mesh, "DG", 0)
	FOS = do.Function(P0)
	FOS.rename("Factor of Safety", "-")

	FOS_vtk = do.File(os.path.join(case_folder, "vtk", "FOS", "FOS.pvd"))

	times, FOS_data = compute_FOS(case_folder)

	print(g.mesh.num_cells())
	print(FOS_data.shape)

	for i, time in enumerate(times):
		FOS.vector()[:] = FOS_data[:,i]
		FOS_vtk << (FOS, time)


if __name__ == '__main__':
	main()