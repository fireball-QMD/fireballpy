import unittest
import json
from os.path import join, split
from numpy.testing import assert_allclose

from ase import Atoms
from ase.io import read
from fireballpy.fdata import load_fdata
from fireballpy import Fireball


class Test(unittest.TestCase):
    TOL = 1e-8
    FDATA = 'CHON'
    FILENAME = 'periodic.xyz'
    DATAFILE = 'periodic.json'
    KPTS = (4,4,4)

    def setUp(self):
        load_fdata(self.FDATA)
        self.assertAllClose = assert_allclose  # Better syntax
        self.atoms = read(join(split(__file__)[0], self.FILENAME))
        with open(self.DATAFILE, 'r', encoding='utf-8') as fp:
            self.data = json.load(fp)
        self.dipoles = ['legacy']
        if self.KPTS == (1,1,1):
            self.dipoles += ['improved']

    def test_lowdin(self):
        assert isinstance(self.atoms, Atoms)
        cm = 'lowdin'
        for dp in self.dipoles:
            data = self.data[f'{dp}::{cm}']
            self.atoms.calc = Fireball(fdata=self.FDATA, kpts=self.KPTS, lazy=False,
                                       charges_method=cm, dipole_method=dp, correction=False)

            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['potential_energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['partial_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_shell_charges(), data['shell_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_eigenvalues(), data['eigenvalues'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['forces'], atol=self.TOL)

    def test_weighted_lowdin(self):
        assert isinstance(self.atoms, Atoms)
        cm = 'weighted_lowdin'
        for dp in self.dipoles:
            data = self.data[f'{dp}::{cm}']
            self.atoms.calc = Fireball(fdata=self.FDATA, kpts=self.KPTS, lazy=False,
                                       charges_method=cm, dipole_method=dp, correction=False)

            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['potential_energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['partial_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_shell_charges(), data['shell_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_eigenvalues(), data['eigenvalues'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['forces'], atol=self.TOL)

    def test_mulliken(self):
        assert isinstance(self.atoms, Atoms)
        cm = 'mulliken'
        for dp in self.dipoles:
            data = self.data[f'{dp}::{cm}']
            self.atoms.calc = Fireball(fdata=self.FDATA, kpts=self.KPTS, lazy=False,
                                       charges_method=cm, dipole_method=dp, correction=False)

            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['potential_energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['partial_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_shell_charges(), data['shell_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_eigenvalues(), data['eigenvalues'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['forces'], atol=self.TOL)

    def test_mulliken_dipole(self):
        assert isinstance(self.atoms, Atoms)
        cm = 'mulliken_dipole'
        for dp in self.dipoles:
            data = self.data[f'{dp}::{cm}']
            self.atoms.calc = Fireball(fdata=self.FDATA, kpts=self.KPTS, lazy=False,
                                       charges_method=cm, dipole_method=dp, correction=False)

            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['potential_energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['partial_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_shell_charges(), data['shell_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_eigenvalues(), data['eigenvalues'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['forces'], atol=self.TOL)

    def test_mulliken_dipole_preserving(self):
        assert isinstance(self.atoms, Atoms)
        cm = 'mulliken_dipole_preserving'
        for dp in self.dipoles:
            data = self.data[f'{dp}::{cm}']
            self.atoms.calc = Fireball(fdata=self.FDATA, kpts=self.KPTS, lazy=False,
                                       charges_method=cm, dipole_method=dp, correction=False)

            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['potential_energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['partial_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_shell_charges(), data['shell_charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.calc.get_eigenvalues(), data['eigenvalues'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['forces'], atol=self.TOL)


if __name__ == '__main__':
    unittest.main()
