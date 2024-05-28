import unittest
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
                                       "save", "molecule_atoms.xyz"))
        self.energy = np.load(os.path.join(FILEPATH,
                                           "save", "molecule_etot.npz"))
        self.charges = np.load(os.path.join(FILEPATH,
                                            "save", "molecule_charges.npz"))
        self.forces = np.load(os.path.join(FILEPATH,
                                           "save", "molecule_forces.npz"))

    def test_lowdin_improved(self):
        self.atoms.calc = Fireball(charges_method='lowdin',
                                   dipole='improved')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_0'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_0'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_0'], atol=self.TOL)

    def test_mulliken_improved(self):
        self.atoms.calc = Fireball(charges_method='mulliken',
                                   dipole='improved')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_1'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_1'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_1'], atol=self.TOL)

    def test_npa_improved(self):
        self.atoms.calc = Fireball(charges_method='npa',
                                   dipole='improved')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_2'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_2'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_2'], atol=self.TOL)

    def test_mulliken_dipole_improved(self):
        self.atoms.calc = Fireball(charges_method='mulliken-dipole',
                                   dipole='improved')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_3'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_3'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_3'], atol=self.TOL)

    def test_mulliken_dipole_preserving_improved(self):
        self.atoms.calc = Fireball(charges_method='mulliken-dipole-preserving',
                                   dipole='improved')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_4'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_4'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_4'], atol=self.TOL)

    def test_lowdin_legacy(self):
        self.atoms.calc = Fireball(charges_method='lowdin',
                                   dipole='legacy')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_5'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_5'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_5'], atol=self.TOL)

    def test_mulliken_legacy(self):
        self.atoms.calc = Fireball(charges_method='mulliken',
                                   dipole='legacy')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_6'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_6'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_6'], atol=self.TOL)

    def test_npa_legacy(self):
        self.atoms.calc = Fireball(charges_method='npa',
                                   dipole='legacy')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_7'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_7'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_7'], atol=self.TOL)

    def test_mulliken_dipole_legacy(self):
        self.atoms.calc = Fireball(charges_method='mulliken-dipole',
                                   dipole='legacy')
        self.atoms.get_potential_energy()
        self.atoms.get_charges()
        self.atoms.get_forces()
        np.testing.assert_allclose(np.array(self.atoms.calc.results['energy']),
                                   self.energy['arr_8'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['shell_charges'],
                                   self.charges['arr_8'], atol=self.TOL)
        np.testing.assert_allclose(self.atoms.calc.results['forces'],
                                   self.forces['arr_8'], atol=self.TOL)

    def test_mulliken_dipole_preserving_legacy(self):
        self.atoms.calc = Fireball(charges_method='mulliken-dipole-preserving',
                                   dipole='legacy')
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
