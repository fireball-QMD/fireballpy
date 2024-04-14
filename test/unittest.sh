#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7

import unittest
import sys
import os
import numpy as np
sys.path.append("..")
from fireballpy import Fireball
from ase.io import read
from ase import Atoms


class Test(unittest.TestCase):

    def test_calcular_ETOT(self):
        atoms = Atoms(symbols=['C', 'H', 'H', 'H', 'H'],
              positions=[( 0.0000,  0.0000,  0.0000),
                         ( 0.6533,  0.6533,  0.6533),
                         (-0.6533,  0.6533,  0.6533),
                         ( 0.6533, -0.6533,  0.6533),
                         ( 0.6533,  0.6533, -0.6533)])

        atoms.calc = Fireball()
        ETOT_cal = atoms.get_potential_energy()
        ETOT_esp = -211.48977661132812
        self.assertAlmostEqual(ETOT_cal, ETOT_esp, delta=1e-12)

        charge_cal = atoms.get_charges()
        charge_esp = np.array( [[1.17513645e+000, 0.0000000],
                                [0.815230489, 0.0000000],
                                [0.876984060, 0.0000000],
                                [0.876984060, 0.0000000],
                                [0.876984060, 0.0000000]])

        #np.testing.assert_allclose(charge_cal, charge_esp, atol=1e-8)

        ftot_cal = atoms.get_forces()
        ftot_esp = np.array([[ 4.63551044,  4.63551044,  4.63551044],
                             [ 2.13428283,  2.13428283,  2.13428283],
                             [-2.53685236, -2.11647058, -2.11647058],
                             [-2.11647058, -2.53685236, -2.11647058],
                             [-2.11647058, -2.11647058, -2.53685236]])

        np.testing.assert_allclose(ftot_cal, ftot_esp, atol=1e-8)


if __name__ == '__main__':
    unittest.main()
