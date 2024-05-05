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

atoms.write('save/test01_atoms.xyz')

for idip in [0,1]:
  for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
    print('-------',C,'-----------------')
    atoms.calc = Fireball(charges=C,
                          igamma   = 1,
                          icluster = 1,
                          idipole  = idip,
                          )
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
  
    print('')
  
    np.savez('save/test01_'+C+'_idipole_'+str(idip)+'.npz', ETOT=np.array(ETOT), charge=charge, force=force )


