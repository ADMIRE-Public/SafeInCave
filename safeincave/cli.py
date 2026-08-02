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

    result = subprocess.run([sys.executable, str(script_path), *args.script_args])
    sys.exit(result.returncode)


def _register_run(subparsers: argparse._SubParsersAction) -> None:
    parser = subparsers.add_parser("run", help="Run <name>.py from the current directory (default: main.py)")
    parser.add_argument(
        "name",
        nargs="?",
        default="main.py",
        help="Script filename, including the .py extension (default: 'main.py')",
    )
    parser.add_argument("script_args", nargs=argparse.REMAINDER, help="Extra arguments forwarded to the script")
    parser.set_defaults(func=_run_script)


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
    _load_plugin_commands(subparsers)

    args = parser.parse_args(sys.argv[1:])
    args.func(args)
