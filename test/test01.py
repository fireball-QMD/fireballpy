import sys
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms


atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
            positions=[( 0.0000, 0.0000, 0.0000),
                       ( 0.6533, 0.6533, 0.6533),
                       (-0.6533, 0.6533, 0.6533),
                       ( 0.6533, -0.6533, 0.6533),
                       ( 0.6533, 0.6533, -0.6533)])


atoms.calc = Fireball(igamma=1,icluster=1)
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