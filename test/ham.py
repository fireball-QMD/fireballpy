import numpy as np
from ase.build import molecule
from fireballpy import Fireball

atoms = molecule('CH4')
atoms.calc = Fireball(fdata='biology')
atoms.get_potential_energy()

norbitals = atoms.calc.norbitals
numorb = np.zeros(len(atoms.calc.atomsystem.species), dtype=np.int64)
hdat = np.zeros((norbitals, norbitals), dtype=np.float64)
sdat = np.zeros((norbitals, norbitals), dtype=np.float64)

from _fireball import get_hs_coords

get_hs_coords(numorb, sdat.T, hdat.T)

print("HAMILTONIAN")
print(hdat - hdat.T)

sdat = 0.5*(sdat + sdat.T)

from scipy.linalg import eigvals

print("\nEIGENVALUES")
hdat = 0.5*(hdat + hdat.T)
eig = np.real(eigvals(hdat, sdat))
eig.sort()
print(eig)

print("\nEIGENVALUES COMPUTED")
print(atoms.calc.get_eigenvalues())
