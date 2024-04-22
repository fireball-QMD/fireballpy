import sys
import os
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms


a=3.569
atoms = Atoms(symbols=['C', 'C'],
            positions=[( 0.0000, 0.0000, 0.0000),
                       (    a/2,    a/2,    a/2)])

atoms.set_cell([(a, a, 0), (a, 0, a), (0, a, a)])

atoms.calc = Fireball(igamma=1, icluster=1, charges="Lowdin") 
#,charges="Mulliken")
# Mulliken-dipole
# NPA:
# Lowdin
# Mulliken-dipole-preserving

ETOT = atoms.get_potential_energy()
print("ETOT = "+str(ETOT))




