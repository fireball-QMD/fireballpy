import numpy as np
from ase.build import molecule
from fireballpy import Fireball

atoms = molecule('CH4')
atoms.calc = Fireball(fdata='biology')
atoms.get_potential_energy()

norbitals = atoms.calc.norbitals
numorb = np.zeros(len(atoms.calc.atomsystem.species), dtype=np.int64)
rowidx = np.zeros(norbitals**2, dtype=np.int64)
colidx = np.zeros(norbitals**2, dtype=np.int64)
hdat = np.zeros(norbitals**2, dtype=np.float64)
sdat = np.zeros(norbitals**2, dtype=np.float64)

from _fireball import get_hs_coords

ndata = get_hs_coords(norbitals, numorb, rowidx, colidx, hdat, sdat)
rowidx = rowidx[:ndata]
colidx = colidx[:ndata]
hdat = hdat[:ndata]
sdat = sdat[:ndata]

from scipy.sparse import coo_matrix
ham = coo_matrix((hdat, (rowidx-1, colidx-1)), shape=(norbitals, norbitals))
over = coo_matrix((sdat, (rowidx-1, colidx-1)), shape=(norbitals, norbitals))

ham = 0.5*(ham + ham.getH()).toarray()
over = 0.5*(over + over.getH()).toarray()

ortoham = np.linalg.inv(over) @ ham
ortoham = 0.5*(ortoham + ortoham.T)

print("get_eigen")
print(atoms.calc.get_eigenvalues())
print("")
print("from ham")
eigenvalues, _ = np.linalg.eig(ortoham)
print(eigenvalues)
