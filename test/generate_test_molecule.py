from ase import Atoms
from ase.io import write
import sys
import numpy as np
sys.path.append("..")
from fireballpy import Fireball


atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
              positions=[(0.0000, 0.0000, 0.0000),
                         (0.6533, 0.6533, 0.6533),
                         (-0.6533, 0.6533, 0.6533),
                         (0.6533, -0.6533, 0.6533),
                         (0.6533, 0.6533, -0.6533)])


ETOT = []
charge = []
force = []

for dipole in ['improved', 'legacy']:
    for C in ['Lowdin', 'Mulliken', 'NPA', 'Mulliken-dipole',
              'Mulliken-dipole-preserving']:
        print('-------', C, '-----------------')
        atoms.calc = Fireball(charges=C, dipole=dipole)
        atoms.get_potential_energy()
        atoms.get_charges()
        atoms.get_forces()
        ie = np.array(atoms.calc.results['energy'])
        ETOT.append(ie)
        print("ETOT = "+str(ie))

        icharge = atoms.calc.results['shell_charges']
        charge.append(icharge)
        print("------atoms.charges-----------")
        for c in icharge:
            print(c)

        iforce = atoms.calc.results['forces']
        force.append(iforce)
        print("------atoms.forces------------")
        for f in iforce:
            print(f)

        print('')

# Descomenta para generar el test otra vez
#write('save/molecule_atoms.xyz', atoms)
#np.savez('save/molecule_etot.npz', *ETOT)
#np.savez('save/molecule_charges.npz', *charge)
#np.savez('save/molecule_forces.npz', *force)
