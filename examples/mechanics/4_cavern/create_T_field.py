import safeincave as sf
from safeincave.Utils import create_field_elems
import numpy as np
import os

def main():
	grid_path = os.path.join("..", "..", "..", "grids", "cavern_overburden_coarse")
	grid = sf.GridHandlerGMSH("geom", grid_path)

	def T_field_fun(x,y,z):
		km = 1000
		dTdZ = 27/km
		T_surface = 20 + 273
		return T_surface - dTdZ*z
	T0_field = create_field_elems(grid, T_field_fun)

	np.savetxt('T.csv', T0_field, delimiter=',', fmt='%f')

	print(T0_field.numpy())

if __name__ == "__main__":
	main()