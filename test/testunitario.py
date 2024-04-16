import sys
import os
import numpy as np
sys.path.append("..")
from fireballpy import Fireball
from ase import Atoms
import unittest

class Test(unittest.TestCase):

  def test_calcular_ETOT(self):
    atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
       positions=[( 0.0000, 0.0000, 0.0000),
             ( 0.6533, 0.6533, 0.6533),
             (-0.6533, 0.6533, 0.6533),
             ( 0.6533, -0.6533, 0.6533),
             ( 0.6533, 0.6533, -0.6533)])

    atoms.calc = Fireball()
    ETOT_cal = atoms.get_potential_energy()
    ETOT_esp = -211.4897918701172
    self.assertAlmostEqual(ETOT_cal, ETOT_esp, delta=1e-12)

    charge_cal = atoms.get_charges()
    charge_esp = np.array( [[1.1751455069,3.3786849976],
                [0.8152207732,0.     ],
                [0.8769829273,0.     ],
                [0.8769829273,0.     ],
                [0.8769829273,0.     ]])

    np.testing.assert_allclose(charge_cal, charge_esp, atol=1e-8)

    ftot_cal = atoms.get_forces()
    ftot_esp = np.array([[ 4.6354694366, 4.6354694366, 4.6354689598],
               [ 2.1343419552, 2.1343419552, 2.1343419552],
               [-2.5368959904,-2.1164574623,-2.1164574623],
               [-2.1164574623,-2.5368959904,-2.1164574623],
               [-2.1164574623,-2.1164574623,-2.5368959904]])

    #arr=ftot_cal
    #print(np.array2string(arr, precision=10, separator=',', suppress_small=True))
    np.testing.assert_allclose(ftot_cal, ftot_esp, atol=1e-8)


if __name__ == '__main__':
  unittest.main()
