import numpy as np
from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.io import write
import sys
sys.path.append("..")
from fireballpy import Fireball
import time

atoms = Atoms(symbols=['Si', 'Si', 'C', 'C', 'Si', 'Si', 'C', 'C'],
              positions=[
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

shell_charges = np.array([
    [1.27983177, 2.14748201, 0.43222374],
    [1.28430411, 2.14317580, 0.42581400],
    [1.24418286, 2.91730000, 0.00000000],
    [1.26634070, 2.91087970, 0.00000000],
    [1.31627088, 2.11941489, 0.38653081],
    [1.30405486, 2.13342448, 0.39705548],
    [1.20839220, 2.93576285, 0.00000000],
    [1.19923677, 2.94832209, 0.00000000]
])

start_time = time.time()

kpoints = monkhorst_pack([4, 4, 4])
C='Mulliken-dipole-preserving'
atoms.calc = Fireball(charges_method=C,kpts=kpoints,shell_charges=shell_charges)
#atoms.calc = Fireball(charges=C)
atoms.get_potential_energy()
print('ETOT =',atoms.calc.results['energy'],'eV')
print('ETOT = -1037.702899 eV saved')
end_time = time.time()
print(f"Tiempo de ejecución: {end_time - start_time} segundos")
# -1037.70289

atoms.get_charges()
print(atoms.calc.results['shell_charges'])
