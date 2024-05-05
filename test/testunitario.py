import sys
import os
import numpy as np
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms
from ase.io import read
import unittest

class Test(unittest.TestCase):
  TOL=1e-4
  def test01(self):
    print('test01')
    atoms = read('save/test01_atoms.xyz')
    for idip in [0,1]:
      for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
        print(C,' idipole = ',idip)
        atoms.calc = Fireball(charges = C,
                          igamma   = 1,
                          icluster = 1,
                          idipole  = idip,
                          )
        valor_esperado = np.load('save/test01_'+C+'_idipole_'+str(idip)+'.npz')
        
        ETOT = atoms.get_potential_energy()
        np.testing.assert_allclose(np.array(ETOT), valor_esperado['ETOT'], atol=self.TOL)
        
        charge = atoms.get_charges()
        np.testing.assert_allclose(charge, valor_esperado['charge'], atol=self.TOL)

        ftot = atoms.get_forces()
        np.testing.assert_allclose(ftot, valor_esperado['force'], atol=self.TOL)


if __name__ == '__main__':
  unittest.main()
