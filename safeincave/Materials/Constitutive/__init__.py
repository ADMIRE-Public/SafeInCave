# Copyright (c) 2026, The SafeInCave Developers
#
# SPDX-License-Identifier: BSD-3-Clause

"""
Constitutive models module for SafeInCave.

This module contains all constitutive element classes for material behavior modeling:
- Elastic: Spring
- Thermoelastic: Thermoelastic
- Non-elastic mechanisms: Viscoelastic, LinearDashpot, DislocationCreep, 
  PressureSolutionCreep, ViscoplasticDesai, MunsonDawsonCreep, ModifiedCamClayViscoplastic

Note: Models are auto-discovered from this folder. New model files are automatically
imported and exported without manual __init__ edits.
"""

import pkgutil
import importlib

__all__ = []

# Dynamically discover and import all model modules
for importer, modname, ispkg in pkgutil.iter_modules(__path__):
    if not modname.startswith('_'):  # Skip private modules
        module = importlib.import_module(f'.{modname}', package=__name__)
        # Assume class name matches module name
        if hasattr(module, modname):
            class_obj = getattr(module, modname)
            globals()[modname] = class_obj
            __all__.append(modname)
        else:
            # Fallback: export all public names from module
            for name in dir(module):
                if not name.startswith('_'):
                    globals()[name] = getattr(module, name)
                    __all__.append(name)
