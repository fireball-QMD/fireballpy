import warnings

import numpy as np

from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

import sys
sys.path.append("/home/dani/fireballpy/build")  # This needs to be thinked
import fireball as fb


class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.

    Out of the box supports molecules containing H, C, N, O and S.

    Parameters
    ----------
    fdata_path : str
        Path to the folder with all FData files.
        More information
        `here <https://fireball-qmd.github.io/fireball.html>`_.

    Notes
    -----
    If you use this code in a publication please cite us. [1]_

    .. [1] Advances and applications in the FIREBALL ab initio tight-binding
       molecular-dynamics formalism.
       James P. Lewis, Pavel Jelínek, José Ortega, Alexander A. Demkov,
       Daniel G. Trabada, Barry Haycock , Hao Wang, Gary Adams,
       John K. Tomfohr , Enrique Abad, Hong Wang, and David A. Drabold.
       Phys. Status Solidi B 248, No. 9, 1989-2007 (2011)
       DOI 10.1002/pssb.201147259
    """

    implemented_properties = ['energy', 'forces', 'charges']

    ignored_changes = ['initial_magmoms']

    def __init__(self, fdata_path: str = None, **kwargs):
        Calculator.__init__(self, **kwargs)
        fb.loadfdata_from_path(fdata_path)

    # Requisite energies
    def _check_compute(self) -> None:
        if 'energy' not in self.results:
            warnings.warn("Energies not computed. Computing energies",
                          UserWarning)
            self._calculate_energies()

    def _calculate_energies(self) -> None:
        fb.call_scf_loop()
        fb.call_getenergy()
        self.energy = fb.get_etot()
        print('ETOT =',self.energy)
        # Save energy
        self.results['energy'] = self.energy
        self.results['free_energy'] = self.energy

    def _calculate_charges(self) -> None:
        self._check_compute()
        self.charges = fb.info_charges(self.natoms)

        # Save charges
        self.results['charges'] = self.charges

    def _calculate_forces(self) -> None:
        self._check_compute()
        fb.call_getforces()
        self.forces = fb.info_forces(self.natoms)

        # Save forces
        self.results['forces'] = self.forces

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes) -> None:

        Calculator.calculate(self, atoms, properties, system_changes)

        # If the atoms change allocate memory
        if 'numbers' in system_changes:
            self.initialize()

        # Check if charges needed
        if 'charges' in properties:
            self._calculate_charges()

        # Check if energy needed
        if 'energy' in properties:
            self._calculate_energies()

        # Check if forces needed
        if 'forces' in properties:
            self._calculate_forces()

    def initialize(self) -> None:
        fb.set_coords(self.atoms.numbers, self.atoms.positions)
        fb.loadlvs_100()
        fb.loadkpts_gamma()
        fb.call_allocate_system()

        self.natoms = len(self.atoms)
        self.charges = np.empty(self.natoms)
        self.forces = np.empty((self.natoms, 3))
