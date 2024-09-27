from __future__ import annotations
from typing import Optional, Any

from numpy.typing import ArrayLike
from scipy.spatial.distance import pdist
from ase import Atoms  # type: ignore
from ase.calculators.abc import GetOutputsMixin  # type: ignore
from ase.calculators.calculator import (Calculator,  # type: ignore
                                        kpts2sizeandoffsets, all_changes)
from ase.dft.kpoints import monkhorst_pack
import numpy as np

from fireballpy._fdata import get_fdata, check_fdata_path, prep_infodat
from fireballpy._errors import raise_fb_error
from fireballpy._options import get_icharge, get_idipole, get_imixer
from fireballpy._dftd3 import get_dftd3

_loaded_fdata = ('', '')

def remove_redundant_kpts(kpts):
    """Function to get rid of negative copies of kpoints"""

    kpts_unique = [kpts[0, :]]
    for k1 in kpts[1:, :]:
        for k2 in kpts_unique:
            if np.allclose(k1, -k2):
                break
        else:
            kpts_unique.append(k1)
    return np.array(kpts_unique, dtype=np.float64)


class Fireball(Calculator, GetOutputsMixin):
    """Python implementation of Fireball DFT code compatible with ASE API.

    The idea of this DFT code is to make use of precomputed integrals of
    the basis functions (called FData) and to implement self consistency
    over the point charges to gain a massive speed improvement.

    This class has been designed with the idea of being used as part of
    `ASE's API <https://wiki.fysik.dtu.dk/ase/>`_.
    As such, it needs to be attached to a valid :class:`ase.Atoms` object.
    
    For a detailed tutorial on how to use this class please look at
    the :ref:`FireballPy User Guide <user_guide>`.

    Parameters
    ----------
    fdata : str
        Name of the FData to be used. See `here <https://fireball-qmd.github.io/fireball.html>`_
        for a table of available FData. If set to ``'custom'`` then
        an ``fdata_path`` pointing to a local FData folder must be provided.
        Please note that the first time an FData is used it needs to download
        all the necessary files.
    lazy : Optional[bool]
        If set to ``True`` (default) it will only load necessary files for the species
        which are involved in the computation. If set to ``False`` it will load
        the FData for all the species available in it. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    kpts : Optional[ArrayLike]
        Three element array-like structure (like list, tuple or numpy.ndarray) with three
        integer elements. A corresponding Monkhorst Pack will be generated with these indices.
    gamma : Optional[bool]
        If ``True`` k-points may be shifted to include the gamma point. By default ``False``.
    verbose : Optional[bool]
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
        By default ``False``.
    charges_method : Optional[str]
        How the autoconsistency in the charges will be performed.
        Options are:
        - NPA (default): ... [NPA]_
        - Mulliken: ... [Mul]_
        - Mulliken-Dipole: ... [MulDip]_
        - Lowdin: ... [Low]_
    dipole : Optional[str]
        Whether to use (``dipole='improved'``) or not (``dipole='legacy'``)
        the improved dipole description.
        By default this is turned on except for periodic systems as
        it is not yet implemented and thus it will be ignored.
    correction : Optional[str]
        If provided it must be 'dftd3'. It will apply London exchange corrections
        with parameters optimised for the selected FData if available.
    mixer_kws : Optional[dict]
        Dictionary with the mixer parameters.
        More information :ref:`here <mixer>`.
    fdata_path : Optional[str]
        Path to the folder with all FData files.
        This parameter is always ignored except when ``fdata='custom'``.

    Notes
    -----
    If you use this code in a publication please cite us [Lewis2011]_.

    References
    ----------
    .. [NPA] Reference for NPA.
           Must add.

    .. [Mul] Reference for Mulliken.
           Must add.

    .. [MulDip] Reference for Mulliken Dipole.
           Must add.

    .. [Low] Reference for Lowdin.
           Must add.

    .. [Lewis2011] James. P. Lewis, et al.,
                   "Advances and applications in the FIREBALL ab initio
                   tight-binding molecular-dynamics formalism",
                   *Phys. Status Solidi B*, 248(9):1989-2007, 2011.
                   `DOI:10.1002/pssb.201147259
                   <https://doi.org/10.1002/pssb.201147259>`_
    """

    implemented_properties = ['energy', 'forces', 'charges']

    ignored_changes = set(['initial_magmoms'])

    default_parameters = {'lazy': True,
                          'kpts': (1, 1, 1),
                          'gamma': False,
                          'verbose': False,
                          'charges_method': 'npa',
                          'dipole_method': 'improved',
                          'correction': None,
                          'mixer': {'method': 'johnson',
                                    'max_iter': np.int64(200),
                                    'mix_order': np.int64(6),
                                    'beta': np.float64(0.1),
                                    'tol': np.float64(1e-8),
                                    'w0': np.float64(0.01)}
                          }

    def __init__(self, fdata: str, *,
                 lazy: bool = True,
                 kpts: ArrayLike = (1,1,1),
                 gamma: bool = False,
                 verbose: bool = False,
                 charges_method: str = 'npa',
                 dipole: str = 'improved',
                 correction: Optional[str] = None,
                 mixer_kws: Optional[dict[str, Any]] = None,
                 fdata_path: Optional[str] = None,
                 **kwargs):

        super().__init__(**kwargs)
        assert self.parameters is not None
        self._scf_computed = False
        self._dftdx = None
        self.nspin = 1
        self.species = None

        # Set parameters
        self.parameters.fdata = fdata
        self.parameters.kpts = kpts
        self.parameters.gamma = bool(gamma)
        self.parameters.lazy = bool(lazy)
        self.parameters.verbose = bool(verbose)
        self.parameters.charges_method = charges_method
        self.parameters.dipole = dipole
        self.parameters.correction = correction
        if mixer_kws is not None:
            self.parameters.mixer['method'] = mixer_kws.get('method', '')
            for prop in ['max_iter', 'mix_order']:
                if prop in mixer_kws:
                    self.parameters.mixer[prop] = np.int64(mixer_kws[prop])
            for prop in ['beta', 'tol', 'w0']:
                if prop in mixer_kws:
                    self.parameters.mixer[prop] = np.float64(mixer_kws[prop])

        # Check FData
        if self.parameters.fdata == 'custom':
            if fdata_path is None:
                raise ValueError("'fdata_path' must be supplied if using custom FData")
            self.parameters.fdata_path = fdata_path

        # Check kpts
        if not hasattr(kpts, '__len__') or kpts.__len__() != 3:  # type: ignore
            raise ValueError("Parameter 'kpts' must be a 3-element "
                                "array-like structure.")

        # Check correction
        if correction is not None:
            if correction == 'dftd3':
                self._dftdx = get_dftd3(fdata)
            else:
                raise ValueError("Parameter 'correction' must None or 'dftd3'")

        # Set Fireball-like options
        self._options = {'dmethod': get_idipole(self.parameters.dipole_method),
                         'qmethod': get_icharge(self.parameters.charges_method),
                         'mixmethod': get_imixer(self.parameters.mixer['method'])}
        self._options.update(self.parameters.mixer)
        del self._options['method']

    def get_eigenvalues(self, kpt=0, spin=0):
        """Get the eigenvalues of the hamiltonian in electronvolts.

        Parameters
        ----------
        kpt : int
            Index of the k-point to obtain the eigenvalues. Check ``Fireball.get_ibz_k_points()``.
            By default takes the first (``kpt = 0``).
        spin : int
            Spin index for the eigenvalues. Currently dummy as no spin polarisation is implemented.

        Returns
        -------
        NDArray[float]
            Array with the eigenvalues for the requested k-point.

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and eigenvalues were not yet computed.
        IndexError
            If the an out-of-bounds k-point (or spin) index is requested.
        """
        return super().get_eigenvalues(kpt, spin)

    def get_fermi_level(self):
        """Get the Fermi level in electronvolts.

        Returns
        -------
        float
            Fermi energy in units of eV.

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and the Fermi level was not yet computed.
        """
        return super().get_fermi_level()

    def get_ibz_k_points(self):
        """Get the k-points used for computation.
        Note that there may be less than imputed as redundant k-points are deleted.

        Returns
        -------
        NDArray[float]
            Array with the coordinates of each k-point in inverse angstroms

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and k-points were not yet computed.
        """
        return super().get_ibz_k_points()

    def get_k_point_weights(self):
        """Get the weights of the k-points used for computation.
        Note that there may be less than imputed as redundant k-points are deleted.

        Returns
        -------
        NDArray[float]
            Array with the weights of each k-point.

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and k-point weights were not yet computed.
        """
        return super().get_k_point_weights()

    def get_number_of_bands(self):
        """Get the number of bands.

        Returns
        -------
        int
            Number of bands

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and the number of bands was not yet computed.
        """
        return super().get_number_of_bands()

    def get_number_of_shells(self):
        """Get the number of shells.

        Returns
        -------
        int
            Number of shells

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and the number of shells was not yet computed.
        """
        return self._get('shell_charges').shape[1]

    def get_shell_charges(self):
        """Get the charges of each shell if computed.

        Returns
        -------
        NDArray[float]
                Array with the charges of each shell.
                A null charge represents that shell is not applicable to the corresponding atom.

        Raises
        ------
        PropertyNotPresent
            If the computation has not been done and shell charges were not yet computed.
        """
        return self._get('shell_charges')

        #    def plot(self, bandpath: Optional[bandpath] = None,
        #             emin: Optional[float] = None,
        #             emax: Optional[float] = None):
        #
        #        xmin = 0
        #        xmax = self.shell_energies.shape[0]
        #
        #        Efermi = get_efermi()
        #        # Ajustamos Efermi
        #        AE = Efermi-self.shell_energies[0][0]
        #        for ek in self.shell_energies:
        #            for e in ek:
        #                if ((Efermi-e) > 0 and abs(Efermi-e) < abs(AE)):
        #                    AE = Efermi-e
        #
        #        Efermi = Efermi-AE
        #
        #        eigen = self.shell_energies-Efermi
        #
        #        if bandpath is None:
        #            X = np.arange(0, xmax, 1)
        #            plt.xticks([])
        #            plt.grid(False)
        #            plt.xlabel('k-points')
        #        else:
        #            X = np.arange(0, xmax, 1)
        #            kpts, Blx, labels = bandpath.get_linear_kpoint_axis()
        #            ticks = Blx*xmax/kpts[-1]
        #            plt.xticks(ticks, labels)
        #            plt.grid(True)
        #
        #        plt.xlim(xmin, xmax)
        #
        #        if emin is None:
        #            emin = np.min(self.shell_energies)-Efermi
        #
        #        if emax is None:
        #            emax = np.max(self.shell_energies)-Efermi
        #
        #        plt.ylim(emin, emax)
        #
        #        for band in range(self.shell_energies.shape[1]):
        #            plt.plot(X, eigen[:, band])
        #
        #        plt.axhline(y=0.0, color='b', linestyle='--', label='y = 50')
        #
        #        plt.ylabel('Energy (eV)')
        #        plt.title('Band Structure')
        #        plt.show()

    # Requisite energies
    def _check_scf(self) -> None:
        if not self._scf_computed or not self.use_cache:
            self._calculate_charges_and_energies()

    def _calculate_charges_and_energies(self) -> None:
        from _fireball import (scf,  # type: ignore
                               get_sizes,
                               get_charges,
                               get_energies,
                               get_eigenvalues)
        assert self.atoms is not None
        assert self.parameters is not None

        # Lauch computation and address errors
        fb_errno = scf(self.parameters.verbose)
        if fb_errno != 0: raise_fb_error(fb_errno)
        self._scf_computed = True

        # Allocate (f2py demands)
        self.nshells, self.nbands = get_sizes()
        self.charges = np.zeros(self.natoms, dtype=np.float64)
        self.shell_charges = np.zeros((self.natoms, self.nshells), dtype=np.float64)
        self.eigenvalues = np.zeros((self.nspin, self.nkpts, self.nbands))

        # Get computation results
        get_charges(self.charges, self.shell_charges.T)
        self.energy, self.fermi_level = get_energies()
        get_eigenvalues(self.eigenvalues[0,:,:].T)

        # Correction
        if self._dftdx is not None:
            dftdx_atoms = self.atoms.copy()
            dftdx_atoms.calc = self._dftdx
            self.energy += dftdx_atoms.get_potential_energy()

        # Store in results
        self.results['charges'] = self.charges
        self.results['shell_charges'] = self.shell_charges
        self.results['energy'] = self.energy
        self.results['fermi_level'] = self.fermi_level
        self.results['eigenvalues'] = self.eigenvalues

    def _calculate_forces(self) -> None:
        from _fireball import calc_forces, get_forces  # type: ignore
        assert self.atoms is not None

        # Energies need to be computed
        self._check_scf()

        # Compute forces and handle errors
        fb_errno = calc_forces()
        if fb_errno != 0: raise_fb_error(fb_errno)

        # Allocate (f2py demands)
        self.forces = np.zeros((self.natoms, 3), dtype=np.float64)

        # Correction
        if self._dftdx is not None:
            dftdx_atoms = self.atoms.copy()
            dftdx_atoms.calc = self._dftdx
            self.forces += dftdx_atoms.get_forces()

        # Get from modules and save
        get_forces(self.forces.T)
        self.results['forces'] = self.forces

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes) -> None:
        from _fireball import update_coords  # type: ignore
        super().calculate(atoms, properties, system_changes)
        assert atoms is not None

        # Check positions
        if np.any(pdist(atoms.get_positions(wrap=True)) < 1e-6):
            raise ValueError("Atom positions are too close! (tolerance 1e-6)")

        # If the atoms change update them in the module
        if 'numbers' in system_changes:
            self._initialize(atoms)
        elif 'positions' in system_changes:
            update_coords(atoms.get_positions().T)

        # Check if SCF loop
        if any(x in properties for x in ['charges', 'energy']):
            self._calculate_charges_and_energies()
        # Check if forces needed
        if 'forces' in properties:
            self._calculate_forces()

    def _initialize_infodat(self, atoms: Atoms) -> None:
        from _fireball import loadfdata_from_path  # type: ignore
        global _loaded_fdata
        assert self.parameters is not None
        if self.parameters.fdata == 'custom':
            check_fdata_path(self.parameters.fdata_path, atoms.get_chemical_symbols())
        else:
            self.parameters.fdata_path = get_fdata(self.parameters.fdata, atoms.get_chemical_symbols())
        if (self.parameters.fdata_path, self.species) != _loaded_fdata:
            self.species = prep_infodat(self.parameters.fdata_path, atoms.get_chemical_symbols(), self.parameters.lazy)
            self.nspecies = len(self.species)
            loadfdata_from_path(self.parameters.fdata_path)
            _loaded_fdata = (self.parameters.fdata_path, self.species)

    def _initialize_options(self, atoms: Atoms) -> None:
        assert self.parameters is not None
        self._options['dmethod'] = np.int64(1)
        self._options['molecule'] = np.int64(1)
        self._options['gonly'] = np.int64(1)
        if not atoms.cell.array.any():
            self.nkpts = 1
            self.ibz_kpoints = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
            self.kpoint_weights = np.array([1.0], dtype=np.float64)
            self._cell = np.diag(2.0*atoms.get_positions().max(axis=0))
            return
        self._options['dmethod'] = np.int64(0)
        self._options['molecule'] = np.int64(0)
        self._cell = atoms.get_cell().array
        if self.parameters.kpts == (1, 1, 1):
            self.nkpts = 1
            self.ibz_kpoints = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
            self.kpoint_weights = np.array([1.0], dtype=np.float64)
            return
        self._options['gonly'] = np.int64(0)
        ksize, koffset = kpts2sizeandoffsets(self.parameters.kpts, gamma=self.parameters.gamma)
        kpts = monkhorst_pack(ksize) + np.array(koffset)
        kpts = 2.0*np.pi*np.dot(kpts, atoms.cell.reciprocal().array.T)
        self.ibz_kpoints = remove_redundant_kpts(kpts)
        self.nkpts = self.ibz_kpoints.shape[0]
        self.kpoint_weights = np.ones(self.nkpts, dtype=np.float64) / self.nkpts

    def _set_options(self, atoms: Atoms) -> None:
        from _fireball import (set_options,  # type: ignore
                               set_coords,
                               set_cell,
                               set_kpoints)
        set_options(**self._options)
        set_coords(atoms.get_atomic_numbers(), atoms.get_positions().T)
        set_cell(self._cell.T)
        set_kpoints(self.ibz_kpoints.T, self.kpoint_weights)

    def _initialize(self, atoms: Atoms) -> None:
        from _fireball import call_allocate_system  # type: ignore
        self.natoms = len(atoms)
        self._initialize_infodat(atoms)
        self._initialize_options(atoms)
        self._set_options(atoms)
        call_allocate_system()

    def _outputmixin_get_results(self) -> dict[str, Any]:
        return self.results
