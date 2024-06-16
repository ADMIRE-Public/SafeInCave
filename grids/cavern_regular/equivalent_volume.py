import numpy as np
import matplotlib.pyplot as plt

def apply_grey_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='0.92')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("0.85")

def trapezoidal_volume(x, y):
	"""
	This function calculates the volume of a solid of revolution (around y=0 axis) based on the trapezoidal rule.
	"""
	volume = 0.0
	area = 0.0
	n = len(x)
	for i in range(1, n):
		R = 0.5*(y[i] + y[i-1])
		A = np.pi*R**2
		d = x[i] - x[i-1]
		area += R*d
		volume += A*d
	return volume, area

def main():
	coord = np.array([
						[0.000000, 205.189718],
						[7.134146, 221.346389],
						[21.951220, 227.809058],
						[47.743902, 235.887393],
						[68.597561, 251.236230],
						[74.634146, 267.392901],
						[72.987805, 285.165239],
						[57.621951, 301.321909],
						[56.524390, 311.823745],
						[47.195122, 331.211750],
						[45.000000, 344.944920],
						[48.292683, 356.254590],
						[47.743902, 380.489596],
						[42.804878, 393.414933],
						[34.573171, 402.301102],
						[19.756098, 412.802938],
						[10.426829, 424.920441],
						[0.000000, 430.000000]
	])
	z_top = max(coord[:,1])
	z_bottom = min(coord[:,1])
	Lz = z_top - z_bottom

	vol_cavern, area_cavern = trapezoidal_volume(coord[:,1], coord[:,0])
	print("A = ", area_cavern)
	print("V = ", vol_cavern)
	print()

	a = 2*np.pi/3
	b = -Lz*np.pi
	c = 0
	d = vol_cavern
	roots = np.roots([a, b, c, d])
	R_eq = abs(roots[1].real)
	D = Lz - 2*R_eq
	vol_eq = 4*np.pi*R_eq**3/3 + D*np.pi*R_eq**2

	print("D = ", D)
	print("R_eq = ", R_eq)
	print("V_eq = ", vol_eq)
	print()


	fig1, ax = plt.subplots(1, 1, figsize=(4, 4))
	fig1.subplots_adjust(top=0.945, bottom=0.115, left=0.160, right=0.945, hspace=0.35, wspace=0.225)

	# Plot original cavern shape
	ax.plot(coord[:,0], coord[:,1], ".-", color="steelblue", linewidth=2.0, label="Original")

	# Plot equivalent cavern shape
	x = np.linspace(0, R_eq, 100000)
	y1 = z_bottom + R_eq
	y2 = z_bottom + R_eq + D

	f1 = lambda x,yc: yc - np.sqrt(R_eq**2 - x**2)
	ax.plot(x, f1(x,y1), "-", color="lightcoral", linewidth=2.0, label="Equivalent")
	f2 = lambda x,yc: yc + np.sqrt(R_eq**2 - x**2)
	ax.plot(x, f2(x,y2), "-", color="lightcoral", linewidth=2.0)
	ax.plot([R_eq, R_eq], [y1, y2], "-", color="lightcoral", linewidth=2.0)
	ax.set_xlabel("x (m)", size=12, fontname="serif")
	ax.set_ylabel("z (m)", size=12, fontname="serif")

	ax.grid(True)
	ax.axis("equal")
	ax.legend(loc=0, shadow=True, fancybox=True)

	z1 = f2(x,y2)
	z2 = f1(x,y1)[::-1]
	z = np.concatenate((z1, z2))[::-1]
	print(z)
	print()	

	x = np.concatenate((x, x[::-1]))
	print(x)
	print()

	coord_eq = np.zeros((len(x), 2))
	coord_eq[:,0] = x
	coord_eq[:,1] = z
	vol_eq, area_eq = trapezoidal_volume(coord_eq[:,1], coord_eq[:,0])
	print("area_eq = ", area_eq)
	print("vol_eq = ", vol_eq)

	apply_grey_theme(fig1, [ax], transparent=True)

	plt.show()


if __name__ == '__main__':
	main()