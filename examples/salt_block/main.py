import os
import sys
sys.path.append(os.path.join("..", "..", "source"))
from Simulator import Simulator
from Utils import *
from dolfin import *
import copy
import time


def main_0():
	start_0 = time.time()

	# Read input file
	input_file = read_json("input_file_0.json")

	# Build simulator
	sim = Simulator(input_file)

	# Run simulation
	sim.run()


	# Print simulation CPU time
	final = time.time()
	print(f"Time: {final - start_0} seconds.")






if __name__ == '__main__':
	main_0()