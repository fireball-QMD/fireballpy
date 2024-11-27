import numpy as np
from ase.build import molecule
from fireballpy import Fireball
from fireballpy.outputs import get_hamiltonian

atoms = molecule('CH4')
atoms.calc = Fireball(fdata='biology')
atoms.get_potential_energy()

ham = get_hamiltonian(atoms=atoms)
