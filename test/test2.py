#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7

import sys
import os
sys.path.append("/home/dani/fireballpy/fireballpy/")
from fireballpy import Fireball

from ase.io import read

atoms = read("/home/dani/fireballpy/test/input.xyz")
atoms.calc = Fireball(fdata_path="/home/dani/Fdata_HC-new/")
ETOT=atoms.get_potential_energy()
print("ETOT = "+str(ETOT))
atoms.get_charges()
atoms.get_potential_energy()
atoms.get_forces()

#print(atoms.energy)
