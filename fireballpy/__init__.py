"""Fireball calcculator for ASE interface

Right now there are three files:

* 'fireball.py' which contains the ``Fireball`` class, subclassed from
ASE's Calculator class
* 'retro.py' which implements utils for older projects which worked with
the pure Fortran version
"""

from ase.calculators.calculator import register_calculator_class  # type: ignore

# Expose API to the user
from fireballpy.fireball import Fireball as Fireball
from fireballpy.version import version as __version__
# from fireballpy.infodat import InfoDat as InfoDat

# Register the calculator as available for ASE
register_calculator_class("fireball", Fireball)
del register_calculator_class

__all__ = [
    '__version__',
    'Fireball',
    #   'InfoDat',
    'retro',
]


def __dir__():
    return __all__
