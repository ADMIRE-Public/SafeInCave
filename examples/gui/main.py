import safeincave as sf
from safeincave.Utils import read_json

def main():
	input_file = read_json("input_file.json")

	sim = sf.Simulator_GUI(input_file)
	sim.run()

if __name__ == '__main__':
	main()