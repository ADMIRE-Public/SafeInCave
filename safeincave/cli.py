"""Command-line entry point for the ``sic`` command.

Built-in subcommands are registered directly below. External packages can
contribute additional subcommands (e.g. ``safeincave_viewer``'s ``view``
command) without owning the ``sic`` console_script themselves, by declaring
an entry point in the ``sic.commands`` group:

    [project.entry-points."sic.commands"]
    view = "some_package.cli:register"

where ``register(subparsers)`` adds its own ``argparse`` subparser, the same
way the built-in ``run`` command does below.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import warnings
from importlib.metadata import entry_points
from pathlib import Path


def _run_script(args: argparse.Namespace) -> None:
    script_path = Path.cwd() / args.name
    if not script_path.is_file():
        print(f"sic run: no such file: {script_path}", file=sys.stderr)
        sys.exit(1)

    if script_path.suffix.lower() in (".yaml", ".yml"):
        # YAML case definition: transpile to Python and run (safeincave.Transpiler)
        command = [sys.executable, "-m", "safeincave.Transpiler", str(script_path), *args.script_args]
    else:
        command = [sys.executable, str(script_path), *args.script_args]
    result = subprocess.run(command)
    sys.exit(result.returncode)


def _register_run(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "run",
        help="Run <name>.py or <name>.yaml from the current directory (default: main.py)",
    )
    parser.add_argument(
        "name",
        nargs="?",
        default="main.py",
        help="Case filename: a .py script or a .yaml case definition (default: 'main.py')",
    )
    parser.add_argument("script_args", nargs=argparse.REMAINDER, help="Extra arguments forwarded to the script")
    parser.set_defaults(func=_run_script)


def _y2p_case(args: argparse.Namespace) -> None:
    case_path = Path.cwd() / args.name
    if not case_path.is_file():
        print(f"sic y2p: no such file: {case_path}", file=sys.stderr)
        sys.exit(1)

    # Options must precede the positional: everything after the case path is
    # captured as script_args by the transpiler's argument parser.
    command = [sys.executable, "-m", "safeincave.Transpiler", "--convert"]
    if args.output:
        command += ["-o", args.output]
    if args.force:
        command += ["--force"]
    command += [str(case_path)]
    result = subprocess.run(command)
    sys.exit(result.returncode)


def _register_y2p(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser(
        "y2p",
        help="Convert a .yaml case definition to an equivalent Python script",
    )
    parser.add_argument("name", help="YAML case filename")
    parser.add_argument("-o", "--output", help="Output .py path (default: <case>.py)")
    parser.add_argument(
        "--force", action="store_true", help="Overwrite an existing non-generated file"
    )
    parser.set_defaults(func=_y2p_case)


def _load_plugin_commands(subparsers: argparse._SubParsersAction) -> None:
    for entry_point in entry_points(group="sic.commands"):
        try:
            register = entry_point.load()
            register(subparsers)
        except Exception as exc:  # noqa: BLE001 - a broken plugin must not break `sic`
            warnings.warn(f"Failed to load sic plugin command '{entry_point.name}': {exc}", RuntimeWarning)


def main() -> None:
    parser = argparse.ArgumentParser(prog="sic", description="SafeInCave command-line tools")
    subparsers = parser.add_subparsers(dest="command", required=True)

    _register_run(subparsers)
    _register_y2p(subparsers)
    _load_plugin_commands(subparsers)

    args = parser.parse_args(sys.argv[1:])
    args.func(args)
