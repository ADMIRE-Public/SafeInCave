import os
import numpy as np
import matplotlib.pyplot as plt
import safeincave as sf
import safeincave.PostProcessingTools as post

def apply_grey_theme(fig, axes, transparent=True, grid_color="0.92", back_color='0.85'):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color=grid_color)
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor(back_color)

def main():
	results_folder = os.path.join("output", "case_0")

	# Plot pressure schedule
	fig, ax = plt.subplots(1, 1, figsize=(5, 3.5))
	fig.subplots_adjust(top=0.970, bottom=0.135, left=0.138, right=0.980, hspace=0.35, wspace=0.225)


	points, time_list, u_field = post.read_node_vector(os.path.join(results_folder, "u", "u.xdmf"))

	point_A = post.find_closest_point([0,0,1], points)
	point_B = post.find_closest_point([0,1,1], points)
	point_C = post.find_closest_point([1,1,1], points)
	point_D = post.find_closest_point([1,0,1], points)

	uz_A = u_field[:,point_A,2]*1000
	uz_B = u_field[:,point_B,2]*1000
	uz_C = u_field[:,point_C,2]*1000
	uz_D = u_field[:,point_D,2]*1000

	time_list /= 60

	ax.plot(time_list, uz_A, ".-", color="#377eb8", label="Point A")
	ax.plot(time_list, uz_B, ".-", color="#ff7f00", label="Point B")
	ax.plot(time_list, uz_C, ".-", color="#4daf4a", label="Point C")
	ax.plot(time_list, uz_D, ".-", color="#f781bf", label="Point D")
	ax.set_xlabel("Time (minutes)", size=12, fontname="serif")
	ax.set_ylabel("Displacement (mm)", size=12, fontname="serif")
	ax.grid(True)
	ax.legend(loc=0, shadow=True, fancybox=True)
	ax.grid(True, color="0.92")
	ax.set_facecolor("0.85")

	plt.show()


if __name__ == '__main__':
	main()