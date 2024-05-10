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
  def test(self):
    test="test02"
    atoms = read('save/'+test+'_atoms.xyz')
    cell =  np.load('save/test02_cell.npz')
    atoms.set_cell(cell['cell'])
    for iga in [0]:
      kpts_read=np.load('save/test02_kpts_igamma'+str(iga)+'.npz')
      kpts=kpts_read['kpts']
      print(kpts)
      for C in ['Lowdin','Mulliken','NPA','Mulliken-dipole','Mulliken-dipole-preserving']:
        print(' igamma = ',iga,C)
        atoms.calc = Fireball(charges = C,
                          igamma   = iga,
                          icluster = 0,
                          idipole  = 0,
                          kpts_monkhorst_pack_ind=kpts
                          )
        valor_esperado = np.load('save/'+test+'_'+C+'_igamma_'+str(iga)+'.npz')
        
        ETOT = atoms.get_potential_energy()
        np.testing.assert_allclose(np.array(ETOT), valor_esperado['ETOT'], atol=self.TOL)
        
        charge = atoms.get_charges()
        np.testing.assert_allclose(charge, valor_esperado['charge'], atol=self.TOL)

        ftot = atoms.get_forces()
        np.testing.assert_allclose(ftot, valor_esperado['force'], atol=self.TOL)

if __name__ == '__main__':
  unittest.main()
