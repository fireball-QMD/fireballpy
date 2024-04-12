"""Fireball calcculator for ASE interface

Right now there are three files:

* 'fireball.py' which contains the ``Fireball`` class, subclassed from
ASE's Calculator class
* 'infodat.py' which contains the ``InfoDat`` class to store data present in
the classical info.dat files
* 'retro.py' which implements utils for older projects which worked with
the pure Fortran version
"""

__all__ = [
    'Fireball',
    'InfoDat'
    'retro'
]

__version__ = '0.0.1'

from ase.calculators.calculator import register_calculator_class  # type: ignore

# Expose API to the user
from fireballpy.fireball import Fireball as Fireball
from fireballpy.infodat import InfoDat as InfoDat

# Register the calculator as available for ASE
register_calculator_class("fireball", Fireball)

# Cleanup
del register_calculator_class
