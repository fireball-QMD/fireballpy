"""Fireball calculator for ASE interface
"""

import importlib as _importlib

# Expose API to the user
## Documentation
from fireballpy.version import version as __version__
from fireballpy.utils import get_fb_home
from fireballpy.fdata import available_fdatas, FDataFiles, load_fdata
from fireballpy.atoms import AtomSystem
from fireballpy.kpoints import KPoints
from fireballpy.fireball import BaseFireball
from fireballpy.ase import Fireball

# Register the calculator as available for ASE
from ase.calculators.calculator import register_calculator_class  # type: ignore
register_calculator_class("fireball", Fireball)
del register_calculator_class

# Define the new properties
from ase.outputs import _defineprop  # type: ignore
_defineprop('nshells', int)
_defineprop('shell_charges', float, shape=('natoms', 'nshells'))
del _defineprop

submodules = [
    'basis',
    'bands',
    'outputs',
]

__all__ = submodules + [  # type: ignore
    '__version__',
    'get_fb_home',
    'available_fdatas',
    'load_fdata',
    'Fireball',
]


def __dir__():
    return __all__


def __getattr__(name):
    if name in submodules:
        return _importlib.import_module(f'fireballpy.{name}')
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(f"Module 'fireballpy' has no attribute '{name}'")
