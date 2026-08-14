# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Local extension mechanism for SafeInCave.

An *extension* is a directory (an "entry") that declares, via a single-line
``TARGET`` file at its root, which ``safeincave.*`` prefix it provides —
then lays its files out directly at their natural relative path underneath,
no mirroring required:

    extensions/
    └── AutoDiffJAX/
        ├── TARGET                      # contents: "../safeincave/AutoDiffJAX"
        ├── __init__.py                 # -> safeincave.AutoDiffJAX
        └── Constitutive/
            └── Mechanism.py            # -> safeincave.AutoDiffJAX.Constitutive.Mechanism

``TARGET``'s content is the relative path this entry would sit at, relative
to itself, if it were merged straight into ``safeincave/``'s parent
directory — so ``../safeincave/X`` from inside an entry lands where ``X``
would really live. A deeper single namespace works the same way
(``../safeincave/Materials/Constitutive``).  Every ``.py`` file under an
entry replaces the public module at that relative path; files that don't
exist publicly are added as new modules (and participate in package
auto-discovery, e.g. constitutive models become ``safeincave.MyModel`` like
built-ins).

``TARGET`` declares *one* prefix per entry — for an extension spanning
several unrelated namespaces, group its files under the highest prefix they
share and lay out the rest as real subfolders beneath that, e.g. an entry
touching ``Simulation.Convergence``, ``Simulation.Simulators``, and
``Simulation.TimeControl`` declares ``../safeincave/Simulation`` and keeps
``Convergence/``, ``Simulators/``, ``TimeControl/`` as real subfolders
underneath it.

Extensions are discovered exclusively from **installed Python packages**
that register a ``"safeincave.extensions"`` entry point resolving
(``.load()``) to a directory path — that directory is then treated as a
*container* of entries, one subfolder per extension, each with its own
``TARGET``, e.g.:

    extensions/
    ├── ConstitutiveModels/TARGET   # "../safeincave/Materials/Constitutive"
    └── AutoDiffJAX/TARGET          # "../safeincave/AutoDiffJAX"

(or the container path can itself be a single entry, with its own
``TARGET`` directly at its root). No other discovery path exists: extension
code is only active while its providing package is ``pip install``ed in the
current environment; uninstalling it (or never installing it) means
``import safeincave`` behaves exactly like a vanilla checkout. Entries from
different packages are applied in sorted-by-name order, first match wins
per module.

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
    merged straight into ``safeincave/``'s parent directory. Declares which
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


def _add_container_or_entry(found: dict[str, tuple[Path, tuple[str, ...]]], folder: Path) -> None:
    """
    Register ``folder`` as a single entry if it has its own ``TARGET``,
    otherwise treat it as a container directory and register each
    subfolder as its own entry -- an entry-point-supplied directory can be
    either shape.
    """
    if not folder.is_dir():
        return
    if (folder / "TARGET").is_file():
        _add_target_entry(found, folder)
        return
    for entry in sorted(folder.iterdir(), key=lambda p: p.name):
        if entry.name.startswith(".") or not entry.is_dir():
            continue  # skip .gitignore, README.md, hidden entries
        _add_target_entry(found, entry)


def _entry_point_extension_dirs() -> list[Path]:
    """
    Resolve every installed ``"safeincave.extensions"`` entry point to a
    directory ``Path``. Each entry point's value, once loaded, must be a
    ``Path`` (or ``str``) pointing at a directory -- typically a small
    module-level constant in the installed package (see the module
    docstring). Tolerates a broken/misconfigured entry point with a
    warning rather than failing all of discovery over it.
    """
    from importlib.metadata import entry_points

    dirs: list[Path] = []
    try:
        eps = sorted(entry_points(group="safeincave.extensions"), key=lambda ep: ep.name)
    except Exception as e:
        warnings.warn(f"SafeInCave: could not query extension entry points: {e}", RuntimeWarning)
        return dirs

    for ep in eps:
        try:
            value = ep.load()
        except Exception as e:
            warnings.warn(
                f"SafeInCave extension entry point '{ep.name}' failed to load: {e}",
                RuntimeWarning,
            )
            continue
        path = Path(value)
        if path.is_dir():
            dirs.append(path)
        else:
            warnings.warn(
                f"SafeInCave extension entry point '{ep.name}' resolved to "
                f"'{path}', which is not a directory; skipping.",
                RuntimeWarning,
            )
    return dirs


def _discover_extensions() -> dict[str, tuple[Path, tuple[str, ...]]]:
    """
    Find extension entries from every installed ``"safeincave.extensions"``
    entry point (first match wins per module, see :func:`_add_target_entry`,
    applied in sorted-by-name order across entry points). Each
    entry-point-resolved directory is either one entry itself (has its own
    ``TARGET``) or a real directory holding several — one subfolder per
    extension, each with its own ``TARGET`` (see
    :func:`_add_container_or_entry`). No other discovery source exists: an
    extension is active only while its providing package is installed.
    """
    found: dict[str, tuple[Path, tuple[str, ...]]] = {}
    for folder in _entry_point_extension_dirs():
        _add_container_or_entry(found, folder)
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

        # Package: merge search paths from EVERY matching entry (not just the
        # first) plus the public directory last, so pkgutil auto-discovery in
        # a package's __init__.py (e.g. Materials/Constitutive's model
        # auto-import) sees modules contributed by all of them — multiple
        # entries are allowed to share one TARGET prefix as long as they don't
        # collide on the same leaf filename (that case is still caught by the
        # multi-match warning above, which fires per-module, not per-package).
        pkg_bases = [b for _, b in matches if b.is_dir()]
        public_pkg_dir = self.public_dir.joinpath(*rel)
        ext_init = next((b / "__init__.py" for b in pkg_bases if (b / "__init__.py").is_file()), None)
        public_init = public_pkg_dir / "__init__.py"
        if ext_init is not None:
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
        search = [str(b) for b in pkg_bases]
        if public_pkg_dir.is_dir():
            search.append(str(public_pkg_dir))
        # Record every contributing entry (what actually fed the merged
        # package), not just the first.
        self.active[fullname] = (
            ", ".join(n for n, b in matches if b.is_dir()),
            str(base) + os.sep,
        )
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
