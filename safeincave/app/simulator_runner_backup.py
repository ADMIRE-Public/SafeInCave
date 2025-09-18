# safeincave/app/simulator_runner.py

import sys
import multiprocessing as mp
from threading import Thread


# ---- module-level worker (picklable) ----
def _sim_worker(conn, jsonfilename):
    # Import inside child to avoid inheriting HDF5 state
    from ..Utils import read_json
    from ..Simulators import Simulator_GUI

    class PipeWriter:
        def __init__(self, send):
            self._send = send
        def write(self, text):
            if text:
                self._send(text)
        def flush(self):  # required by file-like API
            pass

    # redirect stdout/stderr to the pipe
    sys_stdout, sys_stderr = sys.stdout, sys.stderr
    sys.stdout = PipeWriter(conn.send)
    sys.stderr = PipeWriter(conn.send)


    try:
        sim = Simulator_GUI(read_json(jsonfilename))
        sim.run()  # ensure sim closes all writers (use context managers in sim)
    finally:
        # restore stdio and close the pipe
        sys.stdout = sys_stdout
        sys.stderr = sys_stderr
        try:
            conn.close()
        except Exception:
            pass


class SimulatorRunner:
    def __init__(self, output_callback):
        self.output_callback = output_callback
        self.process = None
        self.parent_conn = None
        self.child_conn = None
        self.listener_thread = None
        self.jsonfilename = ""

    def setJsonFile(self, filename1):
        self.jsonfilename = filename1

    def run(self):
        if self.process and self.process.is_alive():
            self.stop()

        mp.set_start_method("spawn", force=True)
        ctx = mp.get_context("spawn")  # critical for HDF5 safety
        # one-way is enough (child -> parent)
        self.parent_conn, self.child_conn = ctx.Pipe(duplex=False)

        # IMPORTANT: target is a module-level function
        self.process = ctx.Process(
            target=_sim_worker,
            args=(self.child_conn, self.jsonfilename),
            daemon=False,
        )
        self.process.start()

        # parent must close its copy of the child end so EOF is delivered properly
        try:
            print("\n 5")
            self.child_conn.close()
        except Exception:
            print("\n 6")
            pass

        def listen_pipe():
            try:
                while True:
                    if self.parent_conn.poll(0.1):
                        print("\n 1")
                        self.write(self.parent_conn.recv())
                    if not self.process.is_alive() and not self.parent_conn.poll():
                        print("\n 2")
                        break
            except (EOFError, OSError):
                print("\n 3")
                self.write("\n[⚠️ Pipe closed. Listener exiting.]\n")
            except Exception as e:
                print("\n 4")
                self.write(f"\n[❌ Error in listener thread: {e}]\n")

        self.listener_thread = Thread(target=listen_pipe, daemon=True)
        self.listener_thread.start()

    def stop(self):
        if self.process and self.process.is_alive():
            # prefer cooperative stop inside sim to avoid corrupting files
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
