from ase.build import molecule
from fireballpy import Fireball
from fireballpy.outputs import get_hamiltonian

atoms = molecule('CH4')
atoms.calc = Fireball('biology', verbose=True)

print("=== COORDINATES ===")
for z, c in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
    print(z, c[0], c[1], c[2])
print("")

print("=== SCF LOOP ===")
atoms.get_potential_energy()
print("")

print("=== ENERGY ===")
print(atoms.get_potential_energy())
print("")

print("=== PARTIAL CHARGES ===")
print(atoms.get_charges())
print("")

print("=== SHELL CHARGES ===")
print(atoms.calc.get_shell_charges())
print("")

print("=== FORCES ===")
print(atoms.get_forces())
print("")

print("=== EIGENVALUES ===")
print(atoms.calc.get_eigenvalues())
print("")

print("=== EIGENVECTORS ===")
print(atoms.calc.get_eigenvectors())
print("")

ham = get_hamiltonian(atoms=atoms)

print("=== HAMILTONIAN ===")
print(ham.hamiltonian)
print("")

print("=== OVERLAP ===")
print(ham.overlap)
print("")
