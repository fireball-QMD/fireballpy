#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7

import sys
import os
sys.path.append("/home/dani/fireballpy/fireballpy/")
from fireballpy import Fireball
from ase.io import read

atoms = read("/home/dani/fireballpy/test/input.xyz")
pos=atoms.get_positions()
print("------atoms.positions--------")
print(pos)

atoms.calc = Fireball(fdata_path="/home/dani/Fdata_HC-new/")
ETOT=atoms.get_potential_energy()
print("ETOT = "+str(ETOT))

charge=atoms.get_charges()
print("------atoms.charges-----------")
for c in charge:
  print(c)

force=atoms.get_forces()
print("------atoms.forces------------")
for f in force:
  print(f)
