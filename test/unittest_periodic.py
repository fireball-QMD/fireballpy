import sys
import os
import numpy as np
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms
from ase.io import read
import unittest

class Test(unittest.TestCase):
  TOL=1e-6
  def test_periodic(self):
    print('periodic')
    valor_esperado_etot = np.load('save/periodic_etot.npz')
    valor_esperado_charges = np.load('save/periodic_charges.npz')
    valor_esperado_forces = np.load('save/periodic_forces.npz')
    atoms = read('save/periodic_atoms.xyz')
    atoms.set_cell = np.load('save/periodic_lvs.npz')
    iload=0
    for opt in ['periodic','periodic_gamma']:
      for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
        print(' options = ',opt,', ',C)
        atoms.calc = Fireball(kpts_monkhorst_pack_ind=[4,4,4], options = opt , charges=C )

        ETOT = atoms.get_potential_energy()
        np.testing.assert_allclose(np.array(ETOT), valor_esperado_etot[f'arr_{iload}'], atol=self.TOL)
        
        charge = atoms.get_charges()
        np.testing.assert_allclose(charge, valor_esperado_charges[f'arr_{iload}'], atol=self.TOL)

        ftot = atoms.get_forces()
        np.testing.assert_allclose(ftot, valor_esperado_forces[f'arr_{iload}'], atol=self.TOL)
        iload=iload+1 
 
if __name__ == '__main__':
  unittest.main()
