import sys
import os
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms



atoms = Atoms(symbols=['Si', 'Si', 'C', 'C', 'Si', 'Si', 'C', 'C', 'Si', 'Si', 'C', 'C', 'Si', 'Si', 'C', 'C', 'H', 'H', 'H', 'H'],
            positions=[ (0.000010,  1.544025,  7.999005),
                        (0.000006, -1.074397,  7.731370),
                        (1.555639, -1.553156,  6.777275),
                        (1.555641,  1.530462,  6.772686),
                        (1.555639,  0.071552,  5.604587),
                        (1.555639,  3.086771,  5.674648),
                        (0.000004, -0.002777,  4.485666),
                        (0.000003,  3.115133,  4.525470),
                        (0.000008, -1.566956,  3.355408),
                        (-0.000004, 1.570878,  3.352878),
                        (1.555639, -1.559492,  2.241648),
                        (1.555633,  1.559680,  2.238680),
                        (1.555634,  0.000000,  1.100000),
                        (1.555634,  3.111270,  1.100000),
                        (0.000000,  0.000000,  0.000000),
                        (0.000000,  3.111270,  0.000000),
                        (0.000000, -0.830000, -0.685000),
                        (0.000000,  0.830000, -0.685000),
                        (0.000000,  2.281270, -0.685000),
                        (0.000000,  3.941269, -0.685000)])



atoms.set_cell([(3.11127, 0.00000, 0.0000), 
                (0.00000, 6.22254, 0.0000), 
                (0.00000, 0.00000, 99.0000)])

for C in ['Lowdin']: #,'Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:

  atoms.calc = Fireball(igamma=0, 
                        icluster=0, 
                        charges=C,
                        idipole = 0,
                        kpts_monkhorst_pack_ind=[16,8,1]) 

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



