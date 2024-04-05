import warnings
from typing import Set, Optional

import numpy as np

from ase import Atoms
from ase.calculators.calculator import (Calculator,all_changes)

import sys
import os
sys.path.append("/home/dani/fireballpy/build")
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
        self._numbers = set()  # Keep track of our atoms
        fb.loadfdata_from_path(fdata_path)

    def _calculate_charges(self, atoms: Atoms) -> None:
        self.charges = fb.info_charges(len(atoms))

        # Save energy
        self.results['charges'] = self.charges

    def _calculate_energies(self, atoms: Atoms) -> None:
        fb.call_getenergy()
        self.energy = fb.info_energy()

        # Save energy
        self.results['energy'] = self.energy
        self.results['free_energy'] = self.energy

    def _calculate_forces(self, atoms: Atoms) -> None:
        # Requisite energies
        if 'energies' not in self.results:
            warnings.warn("Energies not computed. Computing energies",
                          UserWarning)
            self._calculate_energies(atoms)

        fb.call_getforces()
        self.forces = fb.info_forces(len(atoms))

        # Save forces
        self.results['forces'] = self.forces

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes) -> None:

        Calculator.calculate(self, atoms, properties, system_changes)

        # If the atoms change load needed Fdata and allocate memory
        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        # Check if charges needed
        if 'charges' in properties:
            self._calculate_charges(atoms)

        # Check if energy needed
        if 'energy' in properties:
            self._calculate_energies(atoms)

        # Check if forces needed
        if 'forces' in properties:
            self._calculate_forces(atoms)

    def initialize(self, atoms: Atoms) -> None:
        numbers = atoms.get_atomic_numbers()
        if set(numbers) - self._numbers:
            self._numbers = set(numbers)
            #fb.set_coords2(len(atoms), atoms.numbers, atoms.positions)
            fb.set_coords(atoms.numbers, atoms.positions)
            fb.loadlvs_100()
            fb.loadkpts_gamma()
            fb.call_allocate_system()

        fb.call_scf_loop()

        self.charges = np.empty(len(atoms))
        self.forces = np.empty((len(atoms), 3))
