# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Local extension mechanism for SafeInCave.

An *extension* is a directory tree that mirrors the ``safeincave/`` package
layout, made visible through the gitignored ``extensions`` path at the
repository root. The simplest setup is a single symlink straight into a
private extension repository:

    SafeInCave/
    ├── safeincave/                # the public package
    └── extensions -> /path/to/SafeInCave_extensions/safeincave    # symlink

    SafeInCave_extensions/
    └── safeincave/                # mirrors the public package layout
        ├── Simulators.py                      # REPLACES safeincave.Simulators
        └── ConstitutiveModels/
            └── MyModel.py                     # ADDS safeincave.ConstitutiveModels.MyModel

Every ``.py`` file in the extension tree replaces the public module at the
same relative path; files that do not exist publicly are added as new modules
(and participate in package auto-discovery, e.g. constitutive models become
``safeincave.MyModel`` like built-ins).

``extensions`` may equally be a real directory holding *several* extension
entries (subdirectories or symlinks, each a mirror — bare or under a
``safeincave/`` subdir); entries are applied in sorted order, first match
wins per module.

The extension is applied by a :data:`sys.meta_path` finder installed at the
very top of ``safeincave/__init__.py`` — before any submodule import — so
every import of an overlaid module, including imports between ``safeincave``
modules themselves, resolves to the extension file. Extension files are
therefore verbatim drop-in modules: relative imports resolve inside
``safeincave`` and no wrapper code is needed.

Rules and escape hatches:

- ``SAFEINCAVE_NO_EXTENSIONS=1`` disables the mechanism (vanilla SafeInCave —
  use it to check whether a bug is upstream's or an extension's).
- ``SAFEINCAVE_EXTENSIONS_DIR`` may list additional extension folders
  (``os.pathsep``-separated) scanned like ``extensions/`` — useful when
  SafeInCave is not installed from a source checkout.
- :func:`discovered_extensions` and :func:`active_extensions` report what was
  found and which modules were actually shadowed/added — log them with
  simulation output for reproducibility.
- Never monkey-patch: replace whole modules through an extension tree, or
  subclass in your own scripts. Overlaying package ``__init__.py`` files is
  possible but discouraged; prefer leaf modules.
"""

from __future__ import annotations

import os
import sys
import warnings
import importlib.abc
import importlib.util
from pathlib import Path

#: Modules that must never be overlaid (the mechanism itself).
_PROTECTED = {"safeincave.extensions"}

_finder: "_ExtensionFinder | None" = None


def _extension_root(entry: Path) -> Path:
    """Overlay root of an extension entry: its ``safeincave/`` subdir if present."""
    nested = entry / "safeincave"
    return nested if nested.is_dir() else entry


def _has_python_files(path: Path) -> bool:
    """Check if path or any subdirectory contains Python files."""
    try:
        for _ in path.rglob("*.py"):
            return True
    except (OSError, RuntimeError):
        # Handle permission errors or symlink traversal issues
        pass
    return False


def _discover_extensions() -> dict[str, Path]:
    """Find extension trees in ``<repo>/extensions/`` and env-var locations."""
    folders = [Path(__file__).resolve().parent.parent / "extensions"]
    extra = os.environ.get("SAFEINCAVE_EXTENSIONS_DIR", "")
    folders += [Path(p) for p in extra.split(os.pathsep) if p]

    found: dict[str, Path] = {}
    for folder in folders:
        if not folder.is_dir():
            continue
        root = _extension_root(folder)
        if _has_python_files(root):
            # The folder (typically `extensions` as a symlink into a private
            # repo) is itself a single mirror root, not a folder of entries.
            resolved = root.resolve()
            name = (
                resolved.parent.name
                if resolved.name == "safeincave"
                else resolved.name
            )
            found.setdefault(name, root)
            continue
        for entry in sorted(folder.iterdir(), key=lambda p: p.name):
            if entry.name.startswith(".") or not entry.is_dir():
                continue  # skip .gitignore, README.md, hidden entries
            if entry.name in found:
                warnings.warn(
                    f"SafeInCave extension name '{entry.name}' appears in more "
                    f"than one extensions folder; using {found[entry.name]}.",
                    RuntimeWarning,
                )
                continue
            found[entry.name] = _extension_root(entry)
    return found


class _ExtensionFinder(importlib.abc.MetaPathFinder):
    """Meta-path finder that redirects ``safeincave.*`` imports to extensions."""

    def __init__(self, roots: dict[str, Path], public_dir: Path):
        self.roots = roots
        self.public_dir = public_dir
        #: fullname -> (extension name, file/dir path) for every module served.
        self.active: dict[str, tuple[str, str]] = {}

    def find_spec(self, fullname, path=None, target=None):
        if not fullname.startswith("safeincave.") or fullname in _PROTECTED:
            return None
        rel = fullname.split(".")[1:]

        # Apply legacy path aliases for backward compatibility
        # Maps old module paths to their new locations after refactoring
        rel_tuple = tuple(rel)
        if rel_tuple == ("ConstitutiveModels",) or len(rel) > 1 and rel[0] == "ConstitutiveModels":
            # Rewrite safeincave.ConstitutiveModels.* -> safeincave.Materials.Constitutive.*
            rel = ["Materials", "Constitutive"] + rel[1:]

        matches: list[tuple[str, Path]] = []
        for name, root in self.roots.items():
            base = root.joinpath(*rel)
            if base.with_suffix(".py").is_file() or base.is_dir():
                matches.append((name, base))
        if not matches:
            return None
        if len(matches) > 1:
            others = ", ".join(m[0] for m in matches[1:])
            warnings.warn(
                f"Multiple SafeInCave extensions provide '{fullname}'; using "
                f"'{matches[0][0]}', ignoring: {others}",
                RuntimeWarning,
            )
        name, base = matches[0]

        pyfile = base.with_suffix(".py")
        if pyfile.is_file():
            # Plain module: replaced (or added) wholesale by the extension file.
            self.active[fullname] = (name, str(pyfile))
            return importlib.util.spec_from_file_location(fullname, pyfile)

        # Package: merge search paths, extension directory first, so both plain
        # imports and pkgutil auto-discovery prefer/see extension modules.
        public_pkg_dir = self.public_dir.joinpath(*rel)
        ext_init = base / "__init__.py"
        public_init = public_pkg_dir / "__init__.py"
        if ext_init.is_file():
            init = ext_init
        elif public_init.is_file():
            init = public_init
        else:
            warnings.warn(
                f"SafeInCave extension '{name}' provides directory '{base}' for "
                f"'{fullname}' but neither it nor the public package has an "
                "__init__.py; ignored.",
                RuntimeWarning,
            )
            return None
        search = [str(base)]
        if public_pkg_dir.is_dir():
            search.append(str(public_pkg_dir))
        # Record the extension directory (what the extension contributed), not
        # the __init__.py actually used, which may be the public one.
        self.active[fullname] = (name, str(base) + os.sep)
        return importlib.util.spec_from_file_location(
            fullname, init, submodule_search_locations=search
        )


def install() -> None:
    """
    Discover extensions and install the import hook.

    Called once at the very top of ``safeincave/__init__.py``. Idempotent;
    a failure never breaks ``import safeincave``.
    """
    global _finder
    if _finder is not None:
        return
    if os.environ.get("SAFEINCAVE_NO_EXTENSIONS"):
        return
    try:
        roots = _discover_extensions()
    except Exception as exc:
        warnings.warn(
            f"SafeInCave extension discovery failed: {exc}", RuntimeWarning
        )
        return
    if not roots:
        return
    _finder = _ExtensionFinder(roots, Path(__file__).resolve().parent)
    sys.meta_path.insert(0, _finder)


def discovered_extensions() -> dict[str, str]:
    """Return ``{extension name: root directory}`` for every extension found."""
    if _finder is None:
        return {}
    return {name: str(root) for name, root in sorted(_finder.roots.items())}


def active_extensions() -> dict[str, str]:
    """
    Return ``{module name: extension file}`` for every replaced/added module.

    Include this (plus package versions) in run logs for reproducibility.
    """
    if _finder is None:
        return {}
    return {mod: file for mod, (_, file) in sorted(_finder.active.items())}
