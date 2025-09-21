# Minimal, headless runner
def main():
    import argparse, os
    # force headless just in case
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
    os.environ.setdefault("PV_OFFSCREEN", "1")
    os.environ.setdefault("DISPLAY", "")

    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    from ..Utils import read_json
    from ..Simulators import Simulator_GUI
    sim = Simulator_GUI(read_json(args.json))  # implement headless branch if needed
    sim.run()

if __name__ == "__main__":
    main()
