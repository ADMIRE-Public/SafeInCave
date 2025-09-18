import os
import sys
from multiprocessing import Process, Pipe
from threading import Thread
from ..Utils import read_json
from ..Simulators import Simulator_GUI


class SimulatorRunner:
    def __init__(self, output_callback):
        self.output_callback = output_callback
        self.process = None
        self.parent_conn, self.child_conn = None, None
        self.listener_thread = None
        self.jsonfilename = ""

    def setJsonFile(self, filename1):
        self.jsonfilename = filename1

    def run(self):
        # Reset connection and process
        if self.process and self.process.is_alive():
            self.stop()

        self.parent_conn, self.child_conn = Pipe()

        def sim_process(conn, jsonfilename):
            # Redirect stdout/stderr to Pipe
            class PipeWriter:
                def write(self, text):
                    conn.send(text)
                def flush(self):
                    pass

            sys.stdout = PipeWriter()
            sys.stderr = PipeWriter()

            sim = Simulator_GUI(read_json(jsonfilename))
            sim.run()

            # Cleanup
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            conn.close()

        self.process = Process(target=sim_process, args=(self.child_conn, self.jsonfilename))
        self.process.start()

        self.process.terminate()

        def listen_pipe():
            try:
                while True:
                    if self.parent_conn.poll(0.1):  # check if there's data
                        msg = self.parent_conn.recv()
                        self.write(msg)
                    if not self.process.is_alive():
                        break
            except (EOFError, OSError):
                self.write("\n[⚠️ Pipe closed. Listener exiting.]\n")
            except Exception as e:
                self.write(f"\n[❌ Error in listener thread: {e}]\n")

        self.listener_thread = Thread(target=listen_pipe, daemon=True)
        self.listener_thread.start()

    def stop(self):
        if self.process and self.process.is_alive():
            self.process.terminate()
            self.process.join()
            self.write("\n❌ Simulation terminated by user.\n")
        if self.listener_thread and self.listener_thread.is_alive():
            self.listener_thread.join(timeout=1)
        self.process = None
        self.listener_thread = None

    def write(self, text):
        if self.output_callback:
            self.output_callback(text)

    def flush(self):
        pass
