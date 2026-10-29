"""Fireball calculator for ASE interface
"""
import importlib as _importlib

# Expose API to the user
## Documentation
from .version import version as __version__
from .utils import get_fb_home
from .fdata import FData
from .options import ChargesMethod, DipoleMethod, MixerMethod
from . import fdata, atoms, kpoints, fireball
from .ase import Fireball

# Register the calculator as available for ASE
from ase.outputs import _defineprop  # ty: ignore[unresolved-import]
from ase.calculators.calculator import register_calculator_class  # ty: ignore[unresolved-import]
register_calculator_class("fireball", Fireball)
del register_calculator_class

# Define the new properties
_defineprop("nshells", int)
_defineprop("shell_charges", float, shape=("natoms", "nshells"))
del _defineprop

submodules = [
    "basis",
    "bands",
    "outputs",
]

__all__ = submodules + [
    "__version__",
    "get_fb_home",
    "FData",
    "ChargesMethod",
    "DipoleMethod",
    "MixerMethod",
    "Fireball",
    "fdata",
    "atoms",
    "kpoints",
    "fireball"
]

def __dir__():
    return __all__

def __getattr__(name):
    if name in submodules:
        return _importlib.import_module(f"fireballpy.{name}")
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(f"Module 'fireballpy' has no attribute '{name}'")
