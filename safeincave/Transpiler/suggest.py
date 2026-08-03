# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Closest-match suggestions for misspelled YAML keywords."""

from __future__ import annotations

import difflib


def closest(word: str, candidates, n: int = 3) -> list:
    """Candidates closest to ``word``, matched case-insensitively.

    Case matters in the API (``nu`` vs ``Nu``, ``E`` vs ``e``) but a
    capitalization slip is exactly the kind of typo a suggestion should catch,
    so matching ignores case while the returned names keep theirs.
    """
    candidates = list(candidates)
    lowered = {}
    for candidate in candidates:
        lowered.setdefault(candidate.lower(), candidate)
    matches = difflib.get_close_matches(word.lower(), list(lowered), n=n)
    return [lowered[match] for match in matches]
