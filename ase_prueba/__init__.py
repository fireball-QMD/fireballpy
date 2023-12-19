"""Fireball calcculator for ASE interface

Right now there are two files:

* 'fireball.py' which contains the ``Fireball`` class, subclassed from
ASE's Calculator class
* 'retro.py' which implements utils for older projects which worked with
the pure Fortran version
"""

__all__ = [
    'Fireball',
    'retro'
]

__version__ = 'dev0.0.1'

from ase.calculators.calculator import register_calculator_class

from .fireball import Fireball

# Register the calculator as available for ASE
register_calculator_class("fireball", Fireball)