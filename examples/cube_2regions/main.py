import os
import sys
sys.path.append(os.path.join("..", "..", "safeincave"))
from Simulator import Simulator
from Utils import read_json


def main_0():
	# Read input file
	input_file = read_json("input_file.json")

	# Build simulator
	sim = Simulator(input_file)

	# # Run simulation
	# sim.run()






if __name__ == '__main__':
	main_0()