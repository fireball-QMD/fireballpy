import numpy as np
from ase import Atoms
from fireballpy import Fireball
from ase.dft.kpoints import bandpath

a = 5.46
atoms = Atoms(symbols=['Si', 'Si'],
              positions=[(0.0000, 0.0000, 0.0000),
                         (a/4,    a/4,    a/4)])
atoms.set_cell([(0, a/2, a/2), (a/2, 0, a/2), (a/2, a/2, 0)])

charges = np.array([[1.3963999, 2.47211926, 0.13149041],
                    [1.39639524, 2.47208974, 0.13150545]])

path = 'WLGXWK'
npoints = 255
B = bandpath(path=path, cell=atoms.cell, npoints=npoints)
kpoints = B.kpts


atoms.calc = Fireball(charges_method='mdp', kpts=kpoints,
                      shell_charges=charges, ifixcharge=1)
atoms.get_potential_energy()


atoms.calc.plot(bandpath=B, emin=0, emax=20)
