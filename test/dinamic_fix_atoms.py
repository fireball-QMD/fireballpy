import time

from ase import Atoms
from ase.md.verlet import VelocityVerlet
from ase.units import fs
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.constraints import FixAtoms
from ase.dft.kpoints import monkhorst_pack

from fireballpy import Fireball

atoms = Atoms(symbols=['Si', 'Si', 'Si', 'Si', 'Si',
                       'Si', 'Si', 'Si', 'Si', 'Si', 'Si',
                       'Si', 'Si', 'Si', 'Si', 'Si', 'Si',
                       'Si', 'Si', 'Si', 'Si', 'Si', 'Si',
                       'Si', 'Si', 'H', 'Si', 'Si', 'Si',
                       'Si', 'H', 'H', 'Si', 'Si', 'H',
                       'H', 'Si', 'Si', 'Si', 'Si', 'Si',
                       'Si', 'H', 'Si', 'Si', 'Si', 'Si',
                       'H', 'H', 'H', 'H', 'Si', 'Si',
                       'Si', 'Si', 'Si', 'Si', 'H', 'H',
                       'Si'],
              positions=[(0.046368,      0.036419,      5.409027),
                         (-3.308507,      9.755219,      5.347741),
                         (0.057406,     19.482344,      5.371565),
                         (3.423953,     13.639856,      5.313367),
                         (0.059663,     15.593381,      5.367930),
                         (0.054103,      7.807438,      5.381465),
                         (0.056986,     11.696614,      5.332811),
                         (0.054350,      3.916303,      5.367104),
                         (3.425110,      9.760709,      5.379114),
                         (-3.300790,     13.640646,      5.341399),
                         (-3.304463,      5.859035,      5.417022),
                         (3.423106,      5.871419,      5.346053),
                         (5.668948,      9.771793,      4.668642),
                         (-1.054078,      1.966846,      4.658697),
                         (-4.452674,     11.698814,      4.662550),
                         (-4.423384,      7.795270,      4.668974),
                         (-1.063106,      9.753529,      4.650953),
                         (2.298160,     11.694820,      4.645532),
                         (-1.066373,     17.538412,      4.655674),
                         (2.308307,     15.600258,      4.661673),
                         (-1.057404,     13.646487,      4.646575),
                         (-1.076907,      5.865070,      4.674407),
                         (2.299955,      3.930857,      4.653152),
                         (2.300320,      7.817357,      4.662657),
                         (-5.565997,      9.740909,      1.507849),
                         (-5.565997,      9.740909,      0.107849),
                         (-4.443317,     11.685454,      2.301705),
                         (-4.443317,      7.796367,      2.301705),
                         (-2.199999,     15.570909,      1.507849),
                         (-2.199999,      3.910909,      1.507849),
                         (-2.199999,     15.570909,      0.107849),
                         (-2.199999,      3.910909,      0.107849),
                         (-2.197954,     11.685454,      1.507849),
                         (-2.197954,      7.796367,      1.507849),
                         (-2.197954,     11.685454,      0.107849),
                         (-2.197954,      7.796367,      0.107849),
                         (-1.077318,     17.515453,      2.301705),
                         (-1.077318,     13.626367,      2.301705),
                         (-1.077318,      5.855455,      2.301705),
                         (-1.077318,      1.966367,      2.301705),
                         (-1.075271,      9.740909,      2.301705),
                         (1.166001,      9.740909,      1.507849),
                         (1.166001,      9.740909,      0.107849),
                         (1.168045,     17.515453,      1.507849),
                         (1.168045,     13.626367,      1.507849),
                         (1.168045,      5.855455,      1.507849),
                         (1.168045,      1.966367,      1.507849),
                         (1.168045,     17.515453,      0.107849),
                         (1.168045,     13.626367,      0.107849),
                         (1.168045,      5.855455,      0.107849),
                         (1.168045,      1.966367,      0.107849),
                         (2.288682,     11.685454,      2.301705),
                         (2.288682,      7.796367,      2.301705),
                         (2.290728,     15.570909,      2.301705),
                         (2.290728,      3.910909,      2.301705),
                         (4.534042,     11.685454,      1.507849),
                         (4.534042,      7.796367,      1.507849),
                         (4.534042,     11.685454,      0.107849),
                         (4.534042,      7.796367,      0.107849),
                         (5.656729,      9.740909,      2.301705)])

atoms.set_cell([(6.73200000,     11.66000000,      0.00000000),
                (-6.73200000,     11.66000000,      0.00000000),
                (0.00000000,      0.00000000,     99.00000000)])

FRAGMENTS = FixAtoms(indices=[24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
                     39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59])
atoms.set_constraint(FRAGMENTS)

kpoints = monkhorst_pack([4, 4, 4])

start_time = time.time()
atoms.calc = Fireball(charges_method='mdp', kpts=kpoints)
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = VelocityVerlet(atoms, trajectory='dinamic.traj',
                     logfile='dinamic.log', timestep=1*fs)
dyn.run(10)

end_time = time.time()
print(f"Tiempo de ejecución: {end_time - start_time} segundos")
print(dyn.get_number_of_steps())
