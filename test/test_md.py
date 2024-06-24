from ase.optimize import MDMin
from fireballpy import Fireball
from ase import Atoms

atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
              positions=[(0.0000, 0.0000, 0.0000),
                         (-0.6533, -0.6533, -0.6533),
                         (-0.6533, 0.6533, 0.6533),
                         (0.6533, -0.6533, 0.6533),
                         (0.6533, 0.6533, -0.6533)])


atoms.calc = Fireball(charges_method='mdp')

ETOT = atoms.get_potential_energy()
print(ETOT)

opt = MDMin(atoms, trajectory='opt.traj', logfile='opt.log')
opt.run(fmax=0.05)
