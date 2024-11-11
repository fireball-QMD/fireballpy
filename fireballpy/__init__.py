"""Fireball calcculator for ASE interface
"""

# Expose API to the user
from fireballpy.fireball import BaseFireball
from fireballpy.ase import Fireball
from fireballpy.fdata import get_fb_home, available_fdatas, new_fdatafiles, FDataFiles
from fireballpy.atoms import new_atomsystem, AtomSystem
from fireballpy.kpoints import new_kpoints, KPoints
from fireballpy.version import version as __version__

# Register the calculator as available for ASE
from ase.calculators.calculator import register_calculator_class  # type: ignore
register_calculator_class("fireball", Fireball)
del register_calculator_class

# Define the new properties
from ase.outputs import _defineprop  # type: ignore
_defineprop('nshells', int)
_defineprop('shell_charges', float, shape=('natoms', 'nshells'))
del _defineprop

__all__ = [
    '__version__',
    'get_fb_home',
    'available_fdatas',
    'Fireball',
]


def __dir__():
    return __all__
