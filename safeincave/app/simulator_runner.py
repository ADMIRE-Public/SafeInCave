import sys, subprocess, threading

class SimulatorRunner:
    def __init__(self, output_callback):
        self.output_callback = output_callback
        self.proc = None
        self.listener_thread = None
        self.jsonfilename = ""

    def setJsonFile(self, filename1):
        self.jsonfilename = filename1

    def run(self):
        self.stop()  # if already running

        # Start a clean Python interpreter running just the CLI
        self.proc = subprocess.Popen(
            [sys.executable, "-m", "safeincave.app.sim_cli", "--json", self.jsonfilename],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            close_fds=True
        )

        def listen():
            try:
                assert self.proc.stdout is not None
                for line in self.proc.stdout:
                    if self.output_callback:
                        self.output_callback(line)
            finally:
                if self.proc and self.proc.stdout:
                    self.proc.stdout.close()

        self.listener_thread = threading.Thread(target=listen, daemon=True)
        self.listener_thread.start()

    def stop(self):
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait()
            if self.output_callback:
                self.output_callback("\n❌ Simulation terminated by user.\n")
        self.proc = None
        # listener thread will exit when stdout closes
