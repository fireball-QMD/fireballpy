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

print("MOLECULE 1")
for idip in [0,1]:
  for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
    print('-------',C,'-----------------')
    atoms.calc = Fireball(charges=C,
                          igamma   = 1,
                          icluster = 1,
                          idipole  = idip,
                          )
    ETOT = atoms.get_potential_energy()
    print("idip = "+str(idip)+"; ETOT = "+str(ETOT))
  
    charge = atoms.get_charges()
#    print("------atoms.charges-----------")
#    for c in charge:
#      print(c)
  
#    force = atoms.get_forces()
#    print("------atoms.forces------------")
#    for f in force:
#      print(f)
  
    print('')

    # Descomenta para generar el test otra vez
    #np.savez('save/test01_'+C+'_idipole_'+str(idip)+'.npz', ETOT=np.array(ETOT), charge=charge, force=force )

a=3.569
atoms = Atoms(symbols=['C', 'C'],
            positions=[( 0.0000, 0.0000, 0.0000),
                       (    a/4,    a/4,    a/4)])
cell=[(a/2, a/2, 0), (a/2, 0, a/2), (0, a/2, a/2)]
atoms.set_cell(cell)

print("MOLECULE 2")
for iga in [0, 1]: #cuando allocatemos y deallo bien ...
  for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
    print('-------',C,'-----------------')

    kpts=[1,1,1]
    if iga == 1:
      kpts=[1,1,1]
      np.savez('save/test02_kpts_igamma'+str(iga)+'.npz', kpts=kpts )
    else:
      kpts=[4,4,4]
      np.savez('save/test02_kpts_igamma'+str(iga)+'.npz', kpts=kpts )

    atoms.calc = Fireball(igamma  = iga, 
                          icluster= 0, 
                          charges = C,
                          idipole = 0,              
                          kpts_monkhorst_pack_ind=kpts) 

    ETOT = atoms.get_potential_energy()
    print("igamma = "+str(iga)+"; ETOT = "+str(ETOT))


    charge = atoms.get_charges()
 #  print("------atoms.charges-----------")
 #  for c in charge:
 #  print(c)

    force = atoms.get_forces()
 #  print("------atoms.forces------------")
 #  for f in force:
 #    print(f)

#    descomenta para volver a generar el test
#    np.savez('save/test02_'+C+'_igamma_'+str(iga)+'.npz', ETOT=np.array(ETOT), charge=charge, force=force )

