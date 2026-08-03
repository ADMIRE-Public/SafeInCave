# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Load a YAML case file, resolving ``!include`` references.

A case can be split across files by tagging any node with ``!include``::

    material: !include materials/salt.yaml

    stages:
      - !include stages/equilibrium.yaml
      - !include stages/operation.yaml

The referenced file's content replaces the tagged node, so
``materials/salt.yaml`` holds the material mapping itself (``density:``,
``elastic:``, ...). Paths are relative to the file containing the tag.

Includes are resolved into plain Python data before any validation runs, so
the rest of the transpiler neither knows nor cares how many files a case came
from.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from .errors import TranspileError


def load(path) -> tuple[object, list[Path]]:
    """Load ``path``, resolving ``!include`` tags.

    Returns the loaded data and every file that contributed to it, the root
    file first.
    """
    sources: list[Path] = []
    data = _load_file(Path(path).resolve(), included_from=None, stack=(), sources=sources)
    return data, sources


def _load_file(path: Path, included_from: Path | None, stack: tuple, sources: list) -> object:
    if path in stack:
        chain = " -> ".join(str(p) for p in (*stack, path))
        raise TranspileError(f"circular !include: {chain}")

    if not path.is_file():
        if included_from is None:
            raise TranspileError(f"{path}: no such file.")
        raise TranspileError(f"{included_from}: !include target not found: {path}")

    sources.append(path)
    loader_cls = _make_loader(path, (*stack, path), sources)
    try:
        with open(path) as f:
            return yaml.load(f, loader_cls)
    except yaml.YAMLError as exc:
        raise TranspileError(f"{path}: invalid YAML: {exc}") from exc


def _make_loader(path: Path, stack: tuple, sources: list):
    """A SafeLoader that resolves ``!include`` relative to ``path``.

    The loader class is built per file so the including file's directory and
    the include chain are bound to it, with no global loader state.
    """

    class IncludeLoader(yaml.SafeLoader):
        pass

    def construct_include(loader, node):
        if not isinstance(node, yaml.ScalarNode):
            raise TranspileError(
                f"{path}: !include takes a single file path. To include several "
                f"files, tag each entry separately, e.g. "
                f"'stages: [!include a.yaml, !include b.yaml]'."
            )
        target = loader.construct_scalar(node)
        if not target:
            raise TranspileError(f"{path}: !include requires a file path.")
        return _load_file((path.parent / target).resolve(), path, stack, sources)

    IncludeLoader.add_constructor("!include", construct_include)
    return IncludeLoader
