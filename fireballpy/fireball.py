from __future__ import annotations
from typing import Optional
import os
import warnings
import numpy as np
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt  # type: ignore
from ase.calculators.calculator import Calculator, all_changes  # type: ignore
from fireballpy.infodat import InfoDat
from fireballpy.fdata import download_needed, get_default_infodat
from ase.dft.kpoints import bandpath


from ._fireball import (call_scf_loop,  # type: ignore
                        call_getenergy,
                        call_getforces,
                        call_allocate_system,
                        loadfdata_from_path,
                        set_coords,
                        set_iqout,
                        set_cell,
                        set_kpoints,
                        set_igamma,
                        set_idipole,
                        set_icluster,
                        set_ifixcharge,
                        set_shell_atom_charge,
                        get_etot,
                        get_efermi,
                        get_nssh,
                        get_atom_force,
                        get_shell_atom_charge,
                        get_fdata_is_load,
                        get_norbitals_new,
                        get_eigen)

ICHARGE_TABLE = {'lowdin': 1, 'mulliken': 2, 'npa': 3,
                 'mulliken_dipole': 4, 'mulliken_dipole_preserving': 7,
                 'md': 4, 'mdp': 7}
IDIPOLE_TABLE = {'improved': 1, 'legacy': 0}
GAMMA = np.array([[0.0, 0.0, 0.0]])
CELL100 = np.array([[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]])


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
    charges_method : Optional[str]
        How the autoconsistency in the charges will be performed.
        Options are:
        - Mulliken (default): ... [1]_
        - Mulliken-Dipole: ... [2]_
        - Mulliken-Dipole-Preserving: ... [3]_
        - Lowdin: ... [4]_
        - NPA: ... [5]_
        Note: this parameter is case insensitive.
    kpts : Optional[ArrayLike]
        K-points considered in the DFT computation.
        These are expressed in units of the unit cell.
        This parameter is only relevant if the Atoms object
        has a defined cell.
        If the Atoms object has a defined cell and this parameter
        is absent, the Gamma point, kpts=[0.0, 0.0, 0.0],
        will be assumed.
    kpts_units : Optional[str]
        Units in which the k-points (if present) should be interpreted.
        There are two options: unit cell units (default),
        kpts_units = 'unit_cell'; and angstroms, kpts_units = 'angstroms'.
    dipole : Optional[str]
        Whether to use (dipole='improved') or not (dipole='legacy')
        the improved dipole description.
        By default this is turned on except for periodic systems as
        it is not yet implemented and thus will be ignored.

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
                 charges_method: Optional[str] = "mulliken",
                 dipole: Optional[str] = "improved",
                 kpts: Optional[ArrayLike] = None,
                 kpts_units: Optional[str] = "unit_cell",
                 shell_charges: Optional[ArrayLike] = None,
                 **kwargs):

        super().__init__(**kwargs)
        self.fdata_path = fdata_path
        self._check_fdata_path()
        self.charges_method = charges_method.lower()
        self._check_charges_method()
        self.dipole = dipole.lower()
        self._check_dipole()
        self.kpts = kpts
        self.kpts_units = kpts_units
        self._check_kpts()
        self.shell_charges = shell_charges
        self._ifixcharges = int(shell_charges is not None)
        
    @staticmethod
    def clean_kpts(kpts):
        kpts_unique = [kpts[0, :]]
        for k1 in kpts[1:, :]:
            for k2 in kpts_unique:
                if np.sum((k1 + k2)**2) < 1e-8:
                    break
            else:
                kpts_unique.append(k1)
        return np.array(kpts_unique)

    def _check_fdata_path(self):
        if self.fdata_path is None:
            return
        if not os.path.isfile(os.path.join(self.fdata_path, 'info.dat')):
            raise ValueError("info.dat file not found in the specified "
                             f"fdata path ({self.fdata_path}).")

    def _check_charges_method(self):
        if self.charges_method not in ICHARGE_TABLE:
            raise ValueError("Parameter 'charges' must be one of "
                             "'lowdin', 'mulliken', 'npa', mulliken_dipole/md "
                             "or 'mulliken_dipole_preserving/mdp'. "
                             f"Got '{self.charges_method}'.")
        self._icharge = ICHARGE_TABLE[self.charges_method]

    def _check_dipole(self):
        if self.dipole not in IDIPOLE_TABLE:
            raise ValueError("Parameter 'dipole' must be either "
                             "'improved' or 'legacy'. "
                             f"Got '{self.dipole}'.")
        self._idipole = IDIPOLE_TABLE[self.dipole]

    def _check_kpts(self):
        if self.kpts_units not in ['unit_cell', 'angstroms']:
            raise ValueError("Parameter 'kpts_units' must be either "
                             "'unit_cell' or 'angstroms'. "
                             f"Got '{self.kpts_units}'.")
        if self.kpts is None:
            return
        self.kpts = np.array(self.kpts)
        if self.kpts.shape[-1] != 3:
            raise ValueError("K-points must have 3 coordinates. "
                             f"Got {self.kpts.shape[-1]}")
        self.kpts = self.kpts.reshape(-1, 3)
        if self.kpts_units == 'unit_cell' and np.max(np.abs(self.kpts)) > 1.0:
            raise ValueError("K-points is not in unit cell coordinates. "
                             "Perhaps try with kpts_units = 'angstroms'.")

    def get_k_point_weights(self):
        w_k=np.ones(len(self.kpts))/len(self.kpts)
        #print(w_k)
        return w_k

    def get_number_of_spins(self):
        return 1             

    def get_eigenvalues(self,kpt,spin):
        self.spin=spin
        index=0
        for ik in range(len(self.kpts)):
           if(np.linalg.norm(self.kpts[ik] - kpt ) < 0.00001):
               index=ik 
        return (self.shell_energies[index]).copy()

    def get_fermi_level(self):
        return 0.0

    def plot(self, bandpath: Optional[bandpath] = None, 
                   emin: Optional[float] = None,
                   emax: Optional[float] = None ):
        
        xmin=0
        xmax=self.shell_energies.shape[0]
        
        Efermi=get_efermi()
        #Ajustamos Efermi
        AE=Efermi-self.shell_energies[0][0]
        for ek in self.shell_energies:
           for e in ek:
             if ((Efermi-e)>0 and abs(Efermi-e)<abs(AE)): 
               AE=Efermi-e

        Efermi=Efermi-AE

        eigen=self.shell_energies-Efermi

        if bandpath is None:
            X = np.arange(0, xmax , 1)
            plt.xticks([])
            plt.grid(False)
            plt.xlabel('k-points')
        else:
            X = np.arange(0, xmax, 1)
            kpts, Blx, labels = bandpath.get_linear_kpoint_axis()
            ticks = Blx*xmax/kpts[-1]
            plt.xticks(ticks, labels)
            plt.grid(True)

        plt.xlim(xmin,xmax)

        if emin is None:
           emin=np.min(self.shell_energies)-Efermi

        if emax is None:
           emax = np.max(self.shell_energies)-Efermi

        plt.ylim(emin,emax)

        for band in range(self.shell_energies.shape[1]):
            plt.plot(X, eigen[:, band])
        
        plt.axhline(y=0.0, color='b', linestyle='--', label='y = 50')

        plt.ylabel('Energy (eV)')
        plt.title('Band Structure')
        plt.show()


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
        self.results['energy'] = self.energy
        self.results['free_energy'] = self.energy
        if not np.allclose(self._kpts, GAMMA):
            self.shell_energies = np.zeros((self._nkpts,
                                            get_norbitals_new()))
            for imu in range(get_norbitals_new()):
                for ik in range(self._nkpts):
                    self.shell_energies[ik, imu] = get_eigen(imu+1, ik+1)
            self.results['shell_energies'] = self.shell_energies

    def _calculate_charges(self) -> None:
        for iatom in range(self.natoms):
            for issh in range(get_nssh(iatom + 1)):
                self.shell_charges[iatom, issh] = \
                    get_shell_atom_charge(issh + 1, iatom + 1)
        self.charges = np.sum(self.shell_charges, axis=1)
        self.results['charges'] = self.charges
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

        if 'positions' in system_changes:
            set_coords(self.atoms.numbers, self.atoms.positions)

        # Check if charges needed
        if 'charges' in properties:
            self._calculate_charges()

        # Check if energy needed
        if 'energy' in properties:
            self._calculate_energies()

        # Check if forces needed
        if 'forces' in properties:
            self._calculate_forces()

    def _initialize_infodat(self):
        if self.fdata_path is None:
            default_infodat = get_default_infodat()
            try:
                self._infodat = default_infodat.select(self.atoms.numbers)
            except KeyError:
                default_infodat = get_default_infodat(refresh=True)
                self._infodat = default_infodat.select(self.atoms.numbers)
            self.fdata_path = download_needed(self._infodat)
        else:
            self._infodat = InfoDat(os.path.join(self.fdata_path, 'info.dat'))
        if get_fdata_is_load() == 0:
            loadfdata_from_path(self.fdata_path)

    def _initialize_options(self):
        if self.atoms.cell.any():  # periodic and peridic_gamma
            self._idipole = 0
            self._icluster = 0
            self._cell = self.atoms.cell.array
            if self.kpts is None:  # periodic_gamma
                warnings.warn("K-points not provided in periodic system."
                              " Gamma point will be assumed", UserWarning)
                self._igamma = 1
                self._kpts = GAMMA
            elif np.allclose(self.kpts, GAMMA):  # periodic_gamma
                self._igamma = 1
                self._kpts = GAMMA
            else:  # periodic
                if self.kpts_units == 'unit_cell':
                    self._kpts = 2*np.pi*np.dot(self.kpts,
                                                self.atoms.cell.reciprocal().T)
                self._kpts = self.clean_kpts(self._kpts)
                self._igamma = 0
        else:  # molecule and molecule_test
            self._icluster = 1
            self._igamma = 1
            self._kpts = GAMMA
            self._cell = CELL100
        self._nkpts = len(self._kpts)
        if self._ifixcharges == 1:
            self.shell_charges = np.array(self.shell_charges)
            if self.shell_charges.shape != (self.natoms, self._infodat.maxshs):
                raise ValueError("Given shell charges have incorrect shape. "
                                 "Expected natoms x max_shells "
                                 f"{(self.natoms, self._infodat.maxshs)}, "
                                 f"got {self.shell_charges.shape}.")

    def _set_options(self):
        set_coords(self.atoms.numbers, self.atoms.positions)
        set_iqout(self._icharge)
        set_idipole(self._idipole)
        set_icluster(self._icluster)
        set_igamma(self._igamma)
        set_kpoints(self._kpts)
        set_cell(self._cell)
        set_ifixcharge(self._ifixcharges)

    def initialize(self) -> None:
        self.natoms = len(self.atoms)
        self._initialize_infodat()
        self._initialize_options()
        self._set_options()

        call_allocate_system()
        self.forces = np.empty((self.natoms, 3))
        if self._ifixcharges == 0:
            self.charges = np.zeros(self.natoms)
            self.shell_charges = np.zeros((self.natoms, self._infodat.maxshs))
        else:
            set_shell_atom_charge(self.shell_charges)
            self.charges = np.sum(self.shell_charges, axis=1)
            self.results['charges'] = self.charges
            self.results['shell_charges'] = self.shell_charges
