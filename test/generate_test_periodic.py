import sys
import os
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms
import numpy as np   

atoms = Atoms(symbols=['Si', 'Si', 'C', 'C', 'Si', 'Si', 'C', 'C'],
    positions = [
    (2.351869, 2.637258, -1.108009),
    (3.459064, 0.714479, 1.110781),
    (1.939962, 1.123183, -2.202729),
    (3.040182, -0.788390, 0.002861),
    (1.502626, -0.402626, -1.099998),
    (2.602627, -2.307882, 1.100000),
    (0.000000, 0.000000, 0.000000),
    (1.100000, -1.905254, 2.200000)
])

atoms.set_cell([
    (1.100000, -1.905255, -2.200000),
    (2.200000, -3.810511, 4.400000),
    (3.810511, 2.200000, 0.000000)
])

atoms.write('save/periodic_atoms.xyz')
np.savez('save/periodic_lvs.npz', atoms.get_cell())


ETOT=[]
charge=[]
force=[]


for opt in ['periodic','periodic_gamma']:
  for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
    print('-------',C,'-----------------')
    #el kpts_monkhorst_pack_ind lo ignora cuando es periodic_gamm
    atoms.calc = Fireball(kpts_monkhorst_pack_ind=[4,4,4], options = opt , charges=C )
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
#np.savez('save/periodic_etot.npz', *ETOT) 
#np.savez('save/periodic_charges.npz', *charge) 
#np.savez('save/periodic_forces.npz', *force)


