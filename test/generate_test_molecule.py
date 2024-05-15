import sys
import numpy as np
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms


atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
            positions=[( 0.0000, 0.0000, 0.0000),
                       ( 0.6533, 0.6533, 0.6533),
                       (-0.6533, 0.6533, 0.6533),
                       ( 0.6533, -0.6533, 0.6533),
                       ( 0.6533, 0.6533, -0.6533)])

atoms.write('save/molecule_atoms.xyz')

ETOT=[]
charge=[]
force=[]

for options in ['molecule','molecule_test']:
  for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
    print('-------',C,'-----------------')
    atoms.calc = Fireball(charges=C, options='molecule')
    ie=atoms.get_potential_energy()
    ETOT.append(np.array(ie))
    print("ETOT = "+str(ie))
  
    icharge = atoms.get_charges()
    charge.append(icharge)
    print("------atoms.charges-----------")
    for c in icharge:
      print(c)
  
    iforce = atoms.get_forces()
    force.append(iforce)
    print("------atoms.forces------------")
    for f in iforce:
      print(f)
  
    print('')

# Descomenta para generar el test otra vez
#np.savez('save/molecule_etot.npz', *ETOT) 
#np.savez('save/molecule_charges.npz', *charge) 
#np.savez('save/molecule_forces.npz', *force) 


