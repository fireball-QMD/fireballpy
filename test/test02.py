import sys
import os
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms


a=3.569
atoms = Atoms(symbols=['C', 'C'],
            positions=[( 0.0000, 0.0000, 0.0000),
                       (    a/4,    a/4,    a/4)])

atoms.set_cell([(a/2, a/2, 0), (a/2, 0, a/2), (0, a/2, a/2)])

atoms.calc = Fireball(igamma=0, icluster=0, 
                      charges="Mulliken-dipole-preserving",
                      idipole = 0,
                      kpts_monkhorst_pack_ind=[4,4,4]) 

ETOT = atoms.get_potential_energy()
print("ETOT = "+str(ETOT))


charge = atoms.get_charges()
print("------atoms.charges-----------")
for c in charge:
    print(c)

force = atoms.get_forces()
print("------atoms.forces------------")
for f in force:
    print(f)



