#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7

from ase.io import read
import sys
sys.path.append("..")
from fireballpy import Fireball

atoms = read("input.xyz")
pos = atoms.get_positions()
print("------atoms.positions--------")
print(pos)

atoms.calc = Fireball(igamma=1)
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
