"""Fireball calcculator for ASE interface

Right now there are three files:

* 'fireball.py' which contains the ``Fireball`` class, subclassed from
ASE's Calculator class
* 'retro.py' which implements utils for older projects which worked with
the pure Fortran version
"""

from ase.calculators.calculator import register_calculator_class  # type: ignore
from ase.outputs import _defineprop  # type: ignore

# Expose API to the user
from fireballpy.ase import Fireball as Fireball
from fireballpy._fdata import get_fb_home, available_fdatas
from fireballpy.version import version as __version__

# Register the calculator as available for ASE
register_calculator_class("fireball", Fireball)
del register_calculator_class

# Define the new properties
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
