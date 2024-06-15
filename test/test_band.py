import numpy as np
from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.dft.kpoints import bandpath
import sys
sys.path.append("..")
from fireballpy import Fireball

a = 5.46
atoms = Atoms(symbols=['Si', 'Si'],
        positions=[( 0.0000, 0.0000, 0.0000),
                       (    a/4,    a/4,    a/4)])  
atoms.set_cell([(0,a/2, a/2), (a/2, 0, a/2), (a/2, a/2,0)])

charges=np.array([[1.396 , 2.472, 0.1315],
                  [1.396 , 2.472, 0.1315]])

print('Si WLGXWK')
path = 'WLGXWK'
npoints = 80
B = bandpath(path=path, cell=atoms.cell, npoints=npoints)
kpoints=B.kpts

atoms.calc = Fireball(charges_method='mdp', kpts=kpoints, shell_charges=charges)
atoms.get_potential_energy()
atoms.calc.plot(bandpath=B,emax=4)
