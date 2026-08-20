# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""Error type raised for any invalid YAML case definition."""


class TranspileError(Exception):
    """Raised when a YAML case definition cannot be converted to Python.

    The message always names the YAML location (section/keyword) that caused
    the failure and, where applicable, lists the legal alternatives, in the
    spirit of keyword validation in conventional FEM input decks.
    """
