from __future__ import annotations
from typing import Optional
import os
import math
import warnings
import numpy as np
from ase.calculators.calculator import Calculator, all_changes  # type: ignore
from fireballpy.infodat import InfoDat
from fireballpy.fdata import download_needed, get_default_infodat
from ase.dft.kpoints import monkhorst_pack
import matplotlib.pyplot as plt
from ._fireball import (call_scf_loop,  # type: ignore
                        call_getenergy,
                        call_getforces,
                        call_allocate_system,
                        loadfdata_from_path,
                        load_cell_100,
                        set_coords,
                        set_iqout,
                        set_cell,
                        set_kpoints,
                        set_igamma,
                        set_idipole,
                        set_icluster,
                        get_etot,
                        get_nssh,
                        get_atom_force,
                        get_eigen,
                        get_norbitals_new,
                        get_shell_atom_charge,
                        get_fdata_is_load)

ICHARGE_TABLE = {'lowdin': 1, 'mulliken': 2, 'npa': 3,
                 'mulliken-dipole': 4,
                 'mulliken-dipole-preserving': 7}
IDIPOLE_TABLE = {'improved': 1, 'legacy': 0}
GAMMA = np.array([[0.0, 0.0, 0.0]])


class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.

    Out of the box supports molecules containing H, B, C, N, O and Si

    Parameters
    ----------
    fdata_path : Optional[str]
        Path to the folder with all FData files. If None (default)
        it will download needed precomputed FData files.
        More information
        `here <https://fireball-qmd.github.io/fireball.html>`_.
    charges : Optional[str]
        How the autoconsistency in the charges will be performed.
        Options are:
        - Mulliken (default): ... [1]_
        - Mulliken-Dipole: ... [2]_
        - Mulliken-Dipole-Preserving: ... [3]_
        - Lowdin: ... [4]_
        - NPA: ... [5]_
        Note: this parameter is case insensitive.
    kpts : Optional[tuple[list[float]] | list[float]]
        K-points considered in the DFT computation.
        These are expressed in units of the unit cell.
        This parameter is only relevant if the Atoms object
        has a defined cell.
        If the Atoms object has a defined cell and this parameter
        is absent, the Gamma point, kpts=[0.0, 0.0, 0.0],
        will be assumed.
    dipole : Optional[str]
        Whether to use (dipole='improved') or not (dipole='legacy')
        the improved dipole description.
        By default this is turned on except for periodic systems as
        it is not yet implemented.

    Notes
    -----
    If you use this code in a publication please cite us. [6]_

    .. [6] Advances and applications in the FIREBALL ab initio tight-binding
       molecular-dynamics formalism.
       James P. Lewis, Pavel Jelínek, José Ortega, Alexander A. Demkov,
       Daniel G. Trabada, Barry Haycock , Hao Wang, Gary Adams,
       John K. Tomfohr , Enrique Abad, Hong Wang, and David A. Drabold.
       Phys. Status Solidi B 248, No. 9, 1989-2007 (2011)
       DOI 10.1002/pssb.201147259
    """

    implemented_properties = ['energy', 'forces', 'charges']

    ignored_changes = ['initial_magmoms']

    def __init__(self, *, fdata_path: Optional[str] = None,
                 charges: Optional[str] = "mulliken",
                 dipole: Optional[str] = "improved",
                 kpts: Optional[tuple[list[float]] | list[float]] = None,
                 **kwargs):

        super().__init__(**kwargs)
        self._fdata_path = fdata_path
        self.charges = charges.lower()
        self.dipole = dipole.lower()
        self.kpts = kpts
        self._check_input()

    def _check_input(self):
        # Charges
        if self.charges not in ICHARGE_TABLE:
            raise ValueError("Parameter 'charges' must be one of "
                             "'lowdin', 'mulliken', 'npa', "
                             "'mulliken-dipole', 'mulliken-dipole-preserving'"
                             f". Got {self.charges}.")
        self._icharge = ICHARGE_TABLE[self.charges]

        # Dipole
        if self.dipole not in IDIPOLE_TABLE:
            raise ValueError("Parameter 'dipole' must be either "
                             "'improved' or 'legacy'"
                             f". Got {self.dipole}.")
        self._idipole = IDIPOLE_TABLE[self.dipole]

    # Requisite energies
    def _check_compute(self) -> None:
        if 'energy' not in self.results:
            warnings.warn(
                "Energies not computed. Computing energies", UserWarning)
            self._calculate_energies()

    def band_structure(self):
        self.energies = np.zeros((len(self.kpts),get_norbitals_new()))
        for imu in range(get_norbitals_new()):
            for ikpoint in range(len(self.kpts)):
                self.energies[ikpoint,imu] = get_eigen(imu+1,ikpoint+1)

        print(self.energies)

    def plot(self):
        plt.figure(figsize=(8, 6))
        k = np.linspace(0, 1, len(self.kpts))
        for band in range(self.energies.shape[1]):
            #plt.plot(self.kpts, self.energies[:, band])
            plt.scatter(k, self.energies[:, band])
        
        plt.xlabel('k-points')
        plt.ylabel('Energy (eV)')
        plt.title('Band Structure')
        plt.grid(True)
        plt.show()

    def _calculate_energies(self) -> None:
        call_scf_loop()
        call_getenergy()
        self.energy = get_etot()
        self.results['energy'] = self.energy
        self.results['free_energy'] = self.energy

    def _calculate_charges(self) -> None:
        for iatom in range(self.natoms):
            for issh in range(get_nssh(iatom + 1)):
                self.shell_charges[iatom, issh] = \
                    get_shell_atom_charge(issh + 1, iatom + 1)
        self.results['charges'] = np.sum(self.shell_charges, axis=1)
        self.results['shell_charges'] = self.shell_charges

    def _calculate_forces(self) -> None:
        self._check_compute()
        call_getforces()
        for iatom in range(self.natoms):
            for i in range(3):
                self.forces[iatom, i] = get_atom_force(iatom+1, i+1)
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

    def kpts_monkhorst_pack_fit_cell(self):
        kpts_all=np.dot(self.kpts,self.atoms.cell.reciprocal().T)*math.pi*2
        aux=[]
        for i in kpts_all:
            poner=True
            for j in aux:
                if np.linalg.norm(i + np.array(j)) < 0.00001:
                    poner=False    
            if(poner):
                aux.append(i)                    
        return np.array(aux)
        


    def initialize(self) -> None:
        self.natoms = len(self.atoms)

        if self._fdata_path is None:
            default_infodat = get_default_infodat()
            self._infodat = default_infodat.select(self.atoms.numbers)
            self._fdata_path = download_needed(self._infodat)
        else:
            self._infodat = InfoDat(os.path.join(self._fdata_path, "info.dat"))

        set_iqout(self._icharge)
        if get_fdata_is_load() == 0:
            loadfdata_from_path(self._fdata_path)

        set_coords(self.atoms.numbers, self.atoms.positions)
        if self.atoms.cell.any():                   #periodic and peridic_gamma
            set_idipole(0)
            set_icluster(0)
            set_cell(self.atoms.cell)
            if self.kpts is None:                   #periodic_gamma
                warnings.warn("K-points not provided in periodic system."
                              " Gamma point will be assumed", UserWarning)
                set_igamma(1)
                self.kpts = GAMMA
                set_kpoints(self.kpts)
            elif np.allclose(self.kpts, GAMMA):     #periodic_gamma
                set_igamma(1)
                self.kpts = GAMMA
                set_kpoints(self.kpts)
            else:                                   #periodic 
                set_igamma(0)
                set_kpoints(self.kpts_monkhorst_pack_fit_cell())

        else:                                       #molecule and molecule_test
            set_idipole(self._idipole)
            set_icluster(1)
            set_igamma(1)
            load_cell_100()
            set_kpoints(GAMMA)

        call_allocate_system()
        self.charges = np.zeros(self.natoms)
        self.shell_charges = np.zeros((self.natoms, self._infodat.maxshs))
        self.forces = np.empty((self.natoms, 3))
