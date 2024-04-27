import sys
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms


for r in [0.0,0.2,0.3]:
  atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
              positions=[(      r, 0.0000, 0.0000),
                         ( 0.6533, 0.6533, 0.6533),
                         (-0.6533, 0.6533, 0.6533),
                         ( 0.6533, -0.6533, 0.6533),
                         ( 0.6533, 0.6533, -0.6533)])
  atoms.calc = Fireball(igamma=1,
                        icluster=1,
                        charges='Mulliken-dipole')
  ETOT = atoms.get_potential_energy()
  print("ETOT = "+str(ETOT))


