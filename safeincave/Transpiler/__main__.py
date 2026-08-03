# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Process entry point for the YAML transpiler.

Invoked by the ``sic`` CLI as ``python -m safeincave.Transpiler``::

    python -m safeincave.Transpiler case.yaml [args...]     # transpile + run
    python -m safeincave.Transpiler case.yaml --convert [-o out.py] [--force]
"""

from __future__ import annotations

import argparse
import sys

from . import TranspileError, convert, run_yaml


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m safeincave.Transpiler",
        description="Convert a YAML case definition to a Python script and run it.",
    )
    parser.add_argument("case", help="Path to the YAML case file")
    parser.add_argument(
        "--convert", action="store_true",
        help="Only write the generated Python script, do not run it",
    )
    parser.add_argument("-o", "--output", help="Output path for --convert")
    parser.add_argument(
        "--force", action="store_true",
        help="With --convert: overwrite an existing non-generated file",
    )
    parser.add_argument(
        "script_args", nargs=argparse.REMAINDER,
        help="Extra arguments forwarded to the generated script",
    )
    args = parser.parse_args(argv)

    try:
        if args.convert:
            out = convert(args.case, args.output, force=args.force)
            print(f"Wrote {out}")
            return 0
        run_yaml(args.case, args.script_args)
        return 0
    except TranspileError as exc:
        print(f"sic: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
