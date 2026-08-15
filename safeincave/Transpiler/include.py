# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Load a single YAML case file.

Splitting a case across files is done with a plain ``{file: <path>}``
mapping, resolved by ``model.py`` at the specific sections that accept it
(``materials``, ``boundaries``, ``loads``, ``outputs``, legacy
``material``) — not here. This module just reads the one file it's given;
see ``model.py``'s ``_resolve_file_ref`` for how ``{file: ...}`` mappings
get expanded, and why that has to happen after this returns rather than
during YAML parsing.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from .errors import TranspileError


def load(path) -> tuple[object, list[Path]]:
    """Load ``path`` as plain YAML.

    Returns the loaded data and the (single-element) list of files that
    contributed to it -- ``model.py`` extends this list as it resolves any
    ``{file: ...}`` references the data contains.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise TranspileError(f"{path}: no such file.")
    sources: list[Path] = [path]
    try:
        with open(path) as f:
            return yaml.safe_load(f), sources
    except yaml.YAMLError as exc:
        raise TranspileError(f"{path}: invalid YAML: {exc}") from exc
