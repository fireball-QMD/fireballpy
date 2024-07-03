import numpy as np
from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.io import write
from fireballpy import Fireball

a= 4.43000*1.5

atoms = Atoms(symbols=['C', 'Si'],
        positions=[( 0.0000, 0.0000, 0.0000),
                   (    a/4,    a/4,    a/4)])  

atoms.set_cell([(a/2, a/2, 0), (a/2, 0, a/2), (0, a/2, a/2)])

ETOT = []
charge = []
force = []

kpoints = monkhorst_pack([4, 4, 4])
for kpts in [None, kpoints]:
    for C in ['Lowdin', 'Mulliken', 'NPA', 'Mulliken_dipole',
              'Mulliken_dipole_preserving']:
        print('-------', C, '-----------------')
        # el kpts_monkhorst_pack_ind lo ignora cuando es periodic_gamm
        atoms.calc = Fireball(charges_method=C, kpts=kpts)
        atoms.get_potential_energy()
        atoms.get_charges()
        atoms.get_forces()
        ie = atoms.calc.results['energy']
        ETOT.append(np.array(ie))
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
write('save/periodic_atoms.xyz', atoms)
#np.savez('save/periodic_etot.npz', *ETOT)
#np.savez('save/periodic_charges.npz', *charge)
#np.savez('save/periodic_forces.npz', *force)
