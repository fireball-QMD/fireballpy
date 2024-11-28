from ase.build import molecule
from fireballpy import Fireball
from fireballpy.outputs import get_hamiltonian

atoms = molecule('CH4')
atoms.calc = Fireball(fdata='biology')
atoms.get_potential_energy()

ham = get_hamiltonian(atoms=atoms)

print("Full matrix")
print(ham.hamiltonian)
print("\nAtom 1: 1; Atom 2: 1; Orbitals 1: 1, 2; Orbitals 2: 2, 3")
print("(0, 0, 1:3, 2:4)")
print(ham.hamiltonian[0,0,1:3,2:4])
