import unittest
from ase.dft.kpoints import monkhorst_pack
from ase.io import read
import sys
import os
import numpy as np
sys.path.append("..")
from fireballpy import Fireball


class Test(unittest.TestCase):
    TOL = 1e-5

    def setUp(self):
        FILEPATH = os.path.dirname(os.path.realpath(__file__))
        self.atoms = read(os.path.join(FILEPATH,
                                       "save", "periodic_atoms.xyz"))
        self.energy = np.load(os.path.join(FILEPATH,
                                           "save", "periodic_etot.npz"))
        self.charges = np.load(os.path.join(FILEPATH,
                                            "save", "periodic_charges.npz"))
        self.forces = np.load(os.path.join(FILEPATH,
                                           "save", "periodic_forces.npz"))
        self.gamma = monkhorst_pack([1, 1, 1])
        self.monkhorst = monkhorst_pack([4, 4, 4])

    def test_lowdin_gamma(self):
        self.atoms.calc = Fireball(charges_method='lowdin',
                                   kpts=self.gamma)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_0'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_0'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_0'], atol=self.TOL)

    def test_mulliken_gamma(self):
        self.atoms.calc = Fireball(charges_method='mulliken',
                                   kpts=self.gamma)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_1'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_1'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_1'], atol=self.TOL)

    def test_npa_gamma(self):
        self.atoms.calc = Fireball(charges_method='npa',
                                   kpts=self.gamma)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_2'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_2'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_2'], atol=self.TOL)

    def test_mulliken_dipole_gamma(self):
        self.atoms.calc = Fireball(charges_method='mulliken_dipole',
                                   kpts=self.gamma)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_3'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_3'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_3'], atol=self.TOL)

    def test_mulliken_dipole_preserving_gamma(self):
        self.atoms.calc = Fireball(charges_method='mulliken_dipole_preserving',
                                   kpts=self.gamma)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_4'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_4'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_4'], atol=self.TOL)

    def test_lowdin_monkhorst(self):
        self.atoms.calc = Fireball(charges_method='lowdin',
                                   kpts=self.monkhorst)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_5'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_5'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_5'], atol=self.TOL)

    def test_mulliken_monkhorst(self):
        self.atoms.calc = Fireball(charges_method='mulliken',
                                   kpts=self.monkhorst)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_6'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_6'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_6'], atol=self.TOL)

    def test_npa_monkhorst(self):
        self.atoms.calc = Fireball(charges_method='npa',
                                   kpts=self.monkhorst)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_7'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_7'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_7'], atol=self.TOL)

    def test_mulliken_dipole_monkhorst(self):
        self.atoms.calc = Fireball(charges_method='mulliken_dipole',
                                   kpts=self.monkhorst)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_8'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_8'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_8'], atol=self.TOL)

    def test_mulliken_dipole_preserving_monkhorst(self):
        self.atoms.calc = Fireball(charges_method='mulliken_dipole_preserving',
                                   kpts=self.monkhorst)
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_9'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_9'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_9'], atol=self.TOL)


if __name__ == '__main__':
    unittest.main()
