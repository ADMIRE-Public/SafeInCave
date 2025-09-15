import os
import sys
import threading

# sys.path.append(os.path.join("..", "..", "safeincave"))
# from Simulator import Simulator
# from Utils import read_json
from ..Simulators import Simulator_GUI

class SimulatorRunner:
    """Class for running the Simulator and managing the output"""
    jsonfilename = "" # "input_file.json"
    def __init__(self, output_callback):
        """
        Directs the Simulator's output to the output_callback function.
        :param output_callback: A function to display the output in the GUI
        """
        self.output_callback = output_callback

    def setJsonFile(self,filename1):
        self.jsonfilename = filename1

    def run(self):
        """Run the simulation in a separate thread"""
        def simulation_thread():
            input_file = read_json(self.jsonfilename)
            sim = SimulatorGUI(input_file)

            # Redirect output to the `output_callback` function
            sys.stdout = self
            sim.run()
            sys.stdout = sys.__stdout__  # Restore standard output

        threading.Thread(target=simulation_thread, daemon=True).start()

    def write(self, text):
        """Receive output and send it to `output_callback`"""
        if self.output_callback:
            self.output_callback(text)

    def flush(self):
        """For compatibility with `sys.stdout`"""
        pass
