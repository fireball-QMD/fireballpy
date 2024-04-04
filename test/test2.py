from ase.io import read
from fireballpy import Fireball

atoms = read("/home/dani/fireballpy/test/input.xyz")
atoms.calc = Fireball(fdata_path="/home/dani/Fdata_HC-new/")

atoms.get_charges()
atoms.get_potential_energy()
atoms.get_forces()
