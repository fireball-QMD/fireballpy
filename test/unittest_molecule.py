import sys
import os
import numpy as np
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms
from ase.io import read
import unittest

FILEPATH = os.path.dirname(os.path.realpath(__file__))

class Test(unittest.TestCase):
  TOL=1e-6
  def test_molecule(self):
    print('molecule')
    valor_esperado_etot = np.load(FILEPATH+'/save/molecule_etot.npz')
    valor_esperado_charges = np.load(FILEPATH+'/save/molecule_charges.npz')
    valor_esperado_forces = np.load(FILEPATH+'/save/molecule_forces.npz')
    atoms = read(FILEPATH+'/save/molecule_atoms.xyz')
    iload=0
    for options in ['molecule','molecule_test']:
      for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
        print(' options = ',options,', ',C)
        atoms.calc = Fireball(charges=C, options='molecule')
        
        ETOT = atoms.get_potential_energy()
        np.testing.assert_allclose(np.array(ETOT), valor_esperado_etot[f'arr_{iload}'], atol=self.TOL)
        
        charge = atoms.get_charges()
        np.testing.assert_allclose(charge, valor_esperado_charges[f'arr_{iload}'], atol=self.TOL)

        ftot = atoms.get_forces()
        np.testing.assert_allclose(ftot, valor_esperado_forces[f'arr_{iload}'], atol=self.TOL)
        iload=iload+1 

if __name__ == '__main__':
  unittest.main()
