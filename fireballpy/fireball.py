from __future__ import annotations
from typing import Optional

import os
import time  # temporal
import warnings

import numpy as np
from ase.calculators.calculator import Calculator, all_changes  # type: ignore

from fireballpy.infodat import InfoDat
from fireballpy.fdata import download_needed, get_default_infodat

from ._fireball import (call_scf_loop,  # type: ignore
                        call_getenergy,
                        call_getforces,
                        call_allocate_system,
                        loadfdata_from_path,
                        set_coords,
                        set_cell,
                        set_iqout,
                        loadlvs_100,
                        loadkpts_gamma,
                        get_etot,
                        get_nssh,
                        get_atom_force,
                        get_shell_atom_charge,
                        set_igamma,
                        set_icluster,
                        get_igamma,
                        get_iqout,
                        get_icluster)


class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.

    Out of the box supports molecules containing H, B, C, N and O

    Parameters
    ----------
    fdata_path : Optional[str]
        Path to the folder with all FData files. If None (default)
        it will download needed precomputed FData files.
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

    def __init__(self, fdata_path: Optional[str] = None, igamma: int = 1,
                 icluster: int = 1, **kwargs):

        super().__init__(**kwargs)
        self._fdata_path = fdata_path
        set_igamma(igamma)
        set_icluster(icluster)

    # Requisite energies
    def _check_compute(self) -> None:
        if 'energy' not in self.results:
            warnings.warn(
                "Energies not computed. Computing energies", UserWarning)
            self._calculate_energies()

    def _calculate_energies(self) -> None:
        call_scf_loop()
        call_getenergy()
        self.energy = get_etot()
        # Save energy
        self.results['energy'] = self.energy
        self.results['free_energy'] = self.energy

    def _calculate_charges(self) -> None:
        for iatom in range(self.natoms):
            for issh in range(get_nssh(iatom + 1)):
                self.charges[iatom, issh] = \
                    get_shell_atom_charge(issh + 1, iatom + 1)
        # Save charges
        self.results['charges'] = self.charges

    def _calculate_forces(self) -> None:
        self._check_compute()
        call_getforces()
        for iatom in range(self.natoms):
            self.forces[iatom, 0] = get_atom_force(iatom+1, 1)
            self.forces[iatom, 1] = get_atom_force(iatom+1, 2)
            self.forces[iatom, 2] = get_atom_force(iatom+1, 3)
        # Save forces
        self.results['forces'] = self.forces

    def calculate(self, atoms=None,
                  properties=['energy'],
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
        self.natoms = len(self.atoms)

        if self._fdata_path is None:
            default_infodat = get_default_infodat()
            self._infodat = default_infodat.select(self.atoms.numbers)
            self._fdata_path = download_needed(self._infodat, self.natoms)
        else:
            self._infodat = InfoDat(os.path.join(self._fdata_path, "info.dat"))

        loadfdata_from_path(self._fdata_path)
        set_coords(self.atoms.numbers, self.atoms.positions)

        if get_icluster() == 1:
            loadlvs_100()
        else:
            set_cell(self.atoms.cell)


        print("get_igamma ", get_igamma())
        if get_igamma() > 0:
            loadkpts_gamma()

        print("get_iquot",get_iqout())
         
        call_allocate_system()

        self.charges = np.zeros((self.natoms, self._infodat.maxshs))
        self.forces = np.empty((self.natoms, 3))
