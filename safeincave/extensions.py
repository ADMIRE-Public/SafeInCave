# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Local extension mechanism for SafeInCave.

An *extension* is a directory (an "entry") that declares, via a single-line
``TARGET`` file at its root, which ``safeincave.*`` prefix it provides —
then lays its files out directly at their natural relative path underneath,
no mirroring required. Made visible through the gitignored ``extensions``
path at the repository root; the simplest setup is a single symlink into a
private extension repository:

    SafeInCave/
    ├── safeincave/                # the public package
    └── extensions -> /path/to/SafeInCave_extensions/extensions   # symlink

    SafeInCave_extensions/extensions/
    └── AutoDiffJAX/
        ├── TARGET                      # contents: "../safeincave/AutoDiffJAX"
        ├── __init__.py                 # -> safeincave.AutoDiffJAX
        └── Constitutive/
            └── Mechanism.py            # -> safeincave.AutoDiffJAX.Constitutive.Mechanism

``TARGET``'s content is the relative path this entry would sit at, relative
to itself, if it were merged straight into a checkout — mirroring the real
``SafeInCave/extensions -> .../extensions`` symlink topology above, where
``extensions/`` and ``safeincave/`` are siblings, so ``../safeincave/X`` from
inside an entry lands where ``X`` would really live. A deeper single
namespace works the same way (``../safeincave/Materials/Constitutive``).
Every ``.py`` file under an entry replaces the public module at that
relative path; files that don't exist publicly are added as new modules
(and participate in package auto-discovery, e.g. constitutive models become
``safeincave.MyModel`` like built-ins).

``TARGET`` declares *one* prefix per entry — for an extension spanning
several unrelated namespaces, group its files under the highest prefix they
share and lay out the rest as real subfolders beneath that, e.g. an entry
touching ``Simulation.Convergence``, ``Simulation.Simulators``, and
``Simulation.TimeControl`` declares ``../safeincave/Simulation`` and keeps
``Convergence/``, ``Simulators/``, ``TimeControl/`` as real subfolders
underneath it.

``extensions`` may itself be a single entry (has its own ``TARGET``) or a
real directory holding *several* entries — one subfolder per extension,
each with its own ``TARGET``, e.g.:

    extensions/
    ├── ConstitutiveModels/TARGET   # "../safeincave/Materials/Constitutive"
    └── AutoDiffJAX/TARGET          # "../safeincave/AutoDiffJAX"

Entries are applied in sorted order, first match wins per module.

The extension is applied by a :data:`sys.meta_path` finder installed at the
very top of ``safeincave/__init__.py`` — before any submodule import — so
every import of an overlaid module, including imports between ``safeincave``
modules themselves, resolves to the extension file. Extension files are
therefore verbatim drop-in modules: relative imports resolve inside
``safeincave`` (by dotted module name, not physical nesting) and no wrapper
code is needed.

Rules and escape hatches:

- Every entry must have a valid ``TARGET``; one without it is skipped with a
  warning (see :func:`_add_target_entry`) rather than guessed at.
- ``SAFEINCAVE_NO_EXTENSIONS=1`` disables the mechanism (vanilla SafeInCave —
  use it to check whether a bug is upstream's or an extension's).
- ``SAFEINCAVE_EXTENSIONS_DIR`` may list additional extension folders
  (``os.pathsep``-separated) — useful when SafeInCave is not installed from a
  source checkout. Each listed path is always a single entry directly (never
  a container of further entries): list one path per extension.
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


def _read_target_prefix(entry: Path) -> tuple[str, ...] | None:
    """
    Read ``entry / "TARGET"`` if present: a single-line relative path of the
    form ``../safeincave/AutoDiffJAX`` (or a deeper prefix,
    ``../safeincave/Materials/Constitutive``, or bare ``../safeincave`` for
    an entry that mirrors several unrelated namespaces and shares no single
    prefix) — the path this entry would sit at, relative to itself, were it
    merged straight into the public checkout (mirroring the real
    ``SafeInCave/extensions -> .../extensions`` symlink topology, where
    ``extensions/`` and ``safeincave/`` are siblings). Declares which
    ``safeincave.*`` namespace this entry's *own root* (no nested
    ``safeincave/``/prefix folders needed) provides — everything after the
    ``safeincave`` path component becomes the dotted prefix (empty for bare
    ``../safeincave``, in which case the entry's own subfolders must mirror
    real ``safeincave.*`` paths in full, same as today's non-empty-prefix
    entries do below whatever prefix they declare).

    Returns ``None`` if there is no ``TARGET`` file. A present-but-invalid
    file (empty, unreadable, or not of the ``.../safeincave[/<prefix>]``
    shape) warns and is treated as absent, matching the module's general
    "discovery never hard-fails" behavior.
    """
    target_file = entry / "TARGET"
    if not target_file.is_file():
        return None
    try:
        text = target_file.read_text(encoding="utf-8").strip()
    except OSError as exc:
        warnings.warn(f"SafeInCave extension '{entry}': could not read TARGET: {exc}", RuntimeWarning)
        return None
    if not text:
        warnings.warn(f"SafeInCave extension '{entry}': TARGET is empty; ignoring.", RuntimeWarning)
        return None

    parts = [p for p in text.replace("\\", "/").split("/") if p and p != "."]
    while parts and parts[0] == "..":
        parts.pop(0)
    if not parts or parts[0] != "safeincave":
        warnings.warn(
            f"SafeInCave extension '{entry}': TARGET must be a relative path of "
            f"the form '../safeincave[/<Prefix>]' (got {text!r}); ignoring.",
            RuntimeWarning,
        )
        return None
    return tuple(parts[1:])


def _add_target_entry(found: dict[str, tuple[Path, tuple[str, ...]]], entry: Path) -> None:
    """
    Register ``entry`` as an extension if it declares a valid ``TARGET``.

    Every extension entry must have one (see :func:`_read_target_prefix`);
    an entry without one is skipped with a warning rather than guessed at —
    discovery never hard-fails, but it also never silently falls back to
    treating an un-marked folder as an extension.
    """
    prefix = _read_target_prefix(entry)
    if prefix is None:
        warnings.warn(
            f"SafeInCave extension folder '{entry}' has no valid TARGET file; "
            "skipping (every extension entry must declare its safeincave.* "
            "prefix via TARGET).",
            RuntimeWarning,
        )
        return
    name = entry.name
    if name in found:
        warnings.warn(
            f"SafeInCave extension name '{name}' appears in more than one "
            f"extensions folder; using {found[name][0]}.",
            RuntimeWarning,
        )
        return
    found[name] = (entry, prefix)


def _discover_extensions() -> dict[str, tuple[Path, tuple[str, ...]]]:
    """Find extension entries in ``<repo>/extensions/`` and env-var locations.

    The default ``extensions`` path is either one entry itself (has its own
    ``TARGET``) or a real directory holding several — one subfolder per
    extension, each with its own ``TARGET`` — applied in sorted order, first
    match wins per module. Each ``SAFEINCAVE_EXTENSIONS_DIR`` entry
    (``os.pathsep``-separated) is always treated as a single entry directly
    — list one path per extension there rather than pointing at a container
    folder.
    """
    found: dict[str, tuple[Path, tuple[str, ...]]] = {}

    default_folder = Path(__file__).resolve().parent.parent / "extensions"
    if default_folder.is_dir():
        if (default_folder / "TARGET").is_file():
            _add_target_entry(found, default_folder)
        else:
            for entry in sorted(default_folder.iterdir(), key=lambda p: p.name):
                if entry.name.startswith(".") or not entry.is_dir():
                    continue  # skip .gitignore, README.md, hidden entries
                _add_target_entry(found, entry)

    extra = os.environ.get("SAFEINCAVE_EXTENSIONS_DIR", "")
    for p in extra.split(os.pathsep):
        if not p:
            continue
        folder = Path(p)
        if folder.is_dir():
            _add_target_entry(found, folder)

    return found


class _ExtensionFinder(importlib.abc.MetaPathFinder):
    """Meta-path finder that redirects ``safeincave.*`` imports to extensions."""

    def __init__(self, roots: dict[str, tuple[Path, tuple[str, ...]]], public_dir: Path):
        self.roots = roots
        self.public_dir = public_dir
        #: fullname -> (extension name, file/dir path) for every module served.
        self.active: dict[str, tuple[str, str]] = {}

    def find_spec(self, fullname, path=None, target=None):
        if not fullname.startswith("safeincave.") or fullname in _PROTECTED:
            return None
        rel = fullname.split(".")[1:]
        rel_tuple = tuple(rel)

        matches: list[tuple[str, Path]] = []
        for name, (root, prefix) in self.roots.items():
            if rel_tuple[: len(prefix)] != prefix:
                continue  # entry declares a TARGET prefix this import isn't under
            base = root.joinpath(*rel[len(prefix):])
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
    return {name: str(root) for name, (root, _prefix) in sorted(_finder.roots.items())}


def active_extensions() -> dict[str, str]:
    """
    Return ``{module name: extension file}`` for every replaced/added module.

    Include this (plus package versions) in run logs for reproducibility.
    """
    if _finder is None:
        return {}
    return {mod: file for mod, (_, file) in sorted(_finder.active.items())}
