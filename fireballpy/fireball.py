from __future__ import annotations
from typing import Optional

import os
import warnings
import numpy as np
from ase.calculators.calculator import Calculator, all_changes  # type: ignore
from ase.dft.kpoints import monkhorst_pack
from fireballpy.infodat import InfoDat
from fireballpy.fdata import download_needed, get_default_infodat

from ._fireball import (call_scf_loop,  # type: ignore
                        call_getenergy,
                        call_getforces,
                        call_allocate_system,
                        loadfdata_from_path,
                        set_coords,
                        set_iqout,
                        set_cell,
                        load_cell_100,
                        set_kpoints,
                        load_kpoints_gamma,
                        set_igamma,
                        set_idipole,
                        set_icluster,
                        get_igamma,
                        get_etot,
                        get_nssh,
                        get_atom_force,
                        get_shell_atom_charge,
                        get_idipole,
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

    _icharge = {'Lowdin': 1, 'Mulliken': 2, 'NPA': 3,
                'Mulliken-dipole': 4, 'Mulliken-dipole-preserving': 7}

    def __init__(self, fdata_path: Optional[str] = None, igamma: int = 1,
                 icluster: int = 1, charges: str = "Mulliken",
                 idipole: int = 0,
                 kpts_monkhorst_pack_ind: list[int] = [1, 1, 1], **kwargs):

        super().__init__(**kwargs)
        self._fdata_path = fdata_path
        self.charges = charges
        set_igamma(igamma)
        set_icluster(icluster)
        set_idipole(idipole)
        self.kpts_monkhorst_pack_ind = kpts_monkhorst_pack_ind

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
            for i in range(3):
                self.forces[iatom, i] = get_atom_force(iatom+1, i+1)
        # Save forces
        self.results['forces'] = self.forces

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes) -> None:
        super().calculate(atoms, properties, system_changes)

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

    def kpts_monkhorst_pack_fb(self):
        kptsunitary = monkhorst_pack(self.kpts_monkhorst_pack_ind)
        kpts_all = np.dot(kptsunitary, self.atoms.cell.reciprocal())*np.pi*2
        aux = []
        for i in kpts_all:
            poner = True
            for j in aux:
                if np.linalg.norm(i + np.array(j)) < 0.00001:
                    poner = False
            if (poner):
                aux.append(i)
        kpts = np.array(aux)
        set_kpoints(kpts)

    def initialize(self) -> None:
        self.natoms = len(self.atoms)

        if self._fdata_path is None:
            default_infodat = get_default_infodat()
            self._infodat = default_infodat.select(self.atoms.numbers)
            self._fdata_path = download_needed(self._infodat)
        else:
            self._infodat = InfoDat(os.path.join(self._fdata_path, "info.dat"))

        loadfdata_from_path(self._fdata_path)
        set_coords(self.atoms.numbers, self.atoms.positions)

        if get_icluster() > 1:
            load_cell_100()
        else:
            set_cell(self.atoms.cell)

        if get_igamma() == 1:
            load_kpoints_gamma()
        else:
            self.kpts_monkhorst_pack_fb()

        set_iqout(self._icharge[self.charges])
        call_allocate_system()

        self.charges = np.zeros((self.natoms, self._infodat.maxshs))
        self.forces = np.empty((self.natoms, 3))
