from __future__ import annotations
from copy import deepcopy
from typing import Optional, Any

from numpy.typing import ArrayLike, NDArray
from scipy.spatial.distance import pdist
from scipy.linalg import solve, pinv, norm
import numpy as np

from fireballpy._fdata import get_fdata, get_fdata_charge_method, check_fdata_path, prep_infodat
from fireballpy._errors import raise_fb_error
from fireballpy._options import get_icharge, get_idipole, get_imixer

from _fireball import (loadfdata_from_path,  # type: ignore
                       set_coords,
                       set_cell,
                       set_kpoints,
                       set_options,
                       call_allocate_system,
                       scf,
                       get_sizes,
                       get_charges,
                       get_energies,
                       get_eigenvalues,
                       calc_forces,
                       get_forces,
                       update_coords)

_loaded_fdata = ('', '')

DEFAULT_MIXER = {'method': 'johnson',
                 'max_iter': np.int64(200),
                 'mix_order': np.int64(6),
                 'beta': np.float64(0.1),
                 'tol': np.float64(1e-8),
                 'w0': np.float64(0.01)}

class BaseFireball:
    """Base Fireball class to be interfaced with any high-level package such as ASE.

    Parameters
    ----------
    fdata : str
        Name of the FData to be used. See ``available_fdatas()``
        for a table of available FData. If set to ``'custom'`` then
        an ``fdata_path`` pointing to a local FData folder must be provided.
        Please note that the first time an FData is used it needs to download
        all the necessary files.
    species : set[str]
        Set with all the species that will appear in the computation.
    numbers : ArrayLike[int]
        Array with the atomic numbers of each atom.
    positions : ArrayLike[float]
        Array with the coordinates of each atom in angstroms.
    wkpts : ArrayLike[float]
        Array with the weights associated to each k-point.
        Does not need to be normalized.
    kpts : ArrayLike[float]
        Array with the k-points in units of the main cell.
        Redundant k-points will be automatically eliminated.
    lazy : bool
        If set to ``True`` it will only load necessary files for the species
        which are involved in the computation. If set to ``False`` it will load
        the FData for all the species available in it. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    verbose bool
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
    charges_method : str
        How the autoconsistency in the charges will be performed.
        For more information see :ref:`here <charges_method>`
    dipole_method : str
        Whether to use (``dipole='improved'``) or not (``dipole='legacy'``)
        the improved dipole description.
        By default this is turned on except for periodic systems as
        it is not yet implemented and thus it will be ignored.
    cell : Optional[ArrayLike[float]]
        Matrix with the cell coordinates. If not provided non-periodicity
        will be assumed.
    pbc : Optional[ArrayLike[bool]]
        Array with information on which direction periodic boundary conditions apply.
        If not provided non-pbc will be assumed.
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

    .. [Mulliken] Reference for Mulliken.
           Must add.

    .. [MullikenDip] Reference for Mulliken Dipole.
           Must add.

    .. [Lowdin] Reference for Lowdin.
           Must add.

    .. [Lewis2011] James. P. Lewis, et al.,
                   "Advances and applications in the FIREBALL ab initio
                   tight-binding molecular-dynamics formalism",
                   *Phys. Status Solidi B*, 248(9):1989-2007, 2011.
                   `DOI:10.1002/pssb.201147259
                   <https://doi.org/10.1002/pssb.201147259>`_
    """

    def __init__(self, *,
                 fdata: str,
                 species: set[str],
                 numbers: ArrayLike,
                 positions: ArrayLike,
                 wkpts: ArrayLike,
                 kpts: ArrayLike,
                 lazy: bool,
                 verbose: bool,
                 charges_method: str = 'auto',
                 dipole_method: str = 'improved',
                 cell: Optional[ArrayLike] = None,
                 pbc: Optional[ArrayLike] = None,
                 mixer_kws: Optional[dict[str, Any]] = None,
                 fdata_path: Optional[str] = None) -> None:

        # Init variables
        self._scf_computed = False
        self._forces_computed = False
        self.nspin = 1
        pbc = [False, False, False] if pbc is None else pbc
        mixer = deepcopy(DEFAULT_MIXER) if mixer_kws is None else deepcopy(mixer_kws)

        # Save to self and ensure type safety
        self.fdata = fdata
        assert isinstance(self.fdata, str)
        self.species = set(species)
        self.nspecies = len(self.species)
        self.numbers = np.ascontiguousarray(numbers, dtype=np.int64)
        assert len(self.numbers.shape) == 1, "'numbers' must be a 1d array"
        self.natoms = self.numbers.size
        self.positions = np.ascontiguousarray(positions, dtype=np.float64)
        assert self.positions.shape == (self.natoms, 3), "'positions' must be a (natoms, 3) array"
        self.wkpts = np.ascontiguousarray(wkpts, dtype=np.float64)
        assert len(self.wkpts.shape) == 1, "'wkpts' must be a 1d array"
        self.nkpts = self.wkpts.size
        self.kpts = np.ascontiguousarray(kpts, dtype=np.float64)
        assert self.kpts.shape == (self.wkpts, 3), "'kpts' must be a (nkpts, 3) array"
        self.lazy = lazy
        assert isinstance(self.lazy, bool)
        self.verbose = verbose
        assert isinstance(self.verbose, bool)
        self.charges_method = charges_method
        assert isinstance(self.charges_method, str)
        self.dipole_method = dipole_method
        assert isinstance(self.dipole_method, str)
        self.cell = np.ascontiguousarray(cell, dtype=np.float64)
        assert self.cell.shape == (3, 3) or np.isnan(self.cell), "'cell' must be a (3, 3) array"
        self.pbc = np.ascontiguousarray(pbc, dtype=np.bool_)
        assert self.pbc.shape == (3,), "'pbc' must be a 3 element array"
        assert isinstance(mixer, dict)
        assert 'method' in mixer, "'mixer' must have a 'method' key"
        self.mixer = {'mixmethod': mixer['method']}
        for prop in ['max_iter', 'mix_order']:
            self.mixer[prop] = np.int64(mixer.get(prop, DEFAULT_MIXER[prop]))
        for prop in ['beta', 'tol', 'w0']:
            self.mixer[prop] = np.float64(mixer.get(prop, DEFAULT_MIXER[prop]))
        self.fdata_path = fdata_path
        if self.fdata == 'custom':
            assert isinstance(self.fdata_path, str), "If 'fdata=\"custom\"' then fdata_path must be provided"

        # Check if too close
        self._wrap_positions()
        if np.any(pdist(self.wrapped_positions) < 1e-6):
            raise ValueError("Atom positions are too close! (tolerance 1e-6)")

        # Load FData and set coordinates
        self._init_fdata()
        set_coords(self.numbers, self.positions.T)

        # Set all periodic/cell variables
        self._init_cell()
        set_cell(self.cell.T)
        set_kpoints(self.kpts.T, self.wkpts)

        # Set Fireball-like options
        self._options = {'dmethod': np.int64(0) if self._isperiodic else get_idipole(dipole_method),
                         'qmethod': get_icharge(self.charges_method),
                         'molecule': np.int64(not self._isperiodic),
                         'gonly': np.int64(self._isgamma)}
        self._options.update(self.mixer)
        self._options['mixmethod'] = get_imixer(self.mixer['mixmethod'])
        set_options(**self._options)

        # Allocate module
        call_allocate_system()

    def _wrap_positions(self, eps: float = 1e-7) -> None:
        if np.isnan(self.cell).any():
            self.wrapped_positions = self.positions
            return
        shift = np.zeros(3, dtype=np.float64)
        shift[self.pbc] -= eps
        fractional = solve(self.cell.T, self.positions.T).T - shift
        for i, periodic in enumerate(self.pbc):
            if periodic:
                fractional[:, i] = shift[i]
        self.wrapped_positions = np.dot(fractional, self.cell)

    def _init_fdata(self) -> None:
        global _loaded_fdata
        if (self.fdata_path, self.species) != _loaded_fdata:
            if self.fdata == 'custom':
                if self.charges_method == 'auto':
                    raise ValueError("Parameter 'charge_method' cannot be 'auto' with custom FData")
                self.fdata_path = check_fdata_path(self.fdata_path, self.species)
            else:
                self.fdata_path = get_fdata(self.fdata, self.species)
                if self.charges_method == 'auto':
                    self.charges_method = get_fdata_charge_method(self.fdata)
            self.species = prep_infodat(self.fdata_path, self.species, self.lazy)
            self.nspecies = len(self.species)
            loadfdata_from_path(self.fdata_path)
            _loaded_fdata = (self.fdata_path, self.species)

    def _remove_redundant_kpts(self) -> None:
        kpts_unique = [self.kpts[0, :]]
        wkpts_unique = [self.wkpts[0]]
        for k1, w1 in zip(self.kpts[1:, :], self.wkpts[1:]):
            for k2 in kpts_unique:
                if norm(k1 + k2) < 1e-5:
                    break
            else:
                kpts_unique.append(k1)
                wkpts_unique.append(w1)
        self.kpts = np.array(kpts_unique, dtype=np.float64)
        self.wkpts = np.array(wkpts_unique, dtype=np.float64)
        self.wkpts /= self.wkpts.sum()
        self.nkpts = self.wkpts.size

    def _init_cell(self) -> None:
        self._isperiodic = not np.isnan(self.cell).any()
        if not self._isperiodic:
            if self.nkpts > 1 or norm(self.kpts[0, :]) > 1e-5:
                raise ValueError("If there is no cell, then 'kpts' must be the gamma point")
            self._isgamma = True
            self.cell = np.diag(2.0*self.positions.max(axis=0))
            self.wkpts /= self.wkpts.sum()
            return
        self.kpts = 2.0*np.pi*np.dot(self.kpts, pinv(self.cell))
        self._remove_redundant_kpts()
        self._isgamma = self.nkpts == 1 and norm(self.kpts[0, :]) < 1e-5

    def run_scf(self) -> None:
        if not self._scf_computed:
            fb_errno = scf(self.verbose)
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            self.nshells, self.nbands = get_sizes()
            self._scf_computed = True

    def calc_forces(self) -> None:
        self.run_scf()
        if not self._forces_computed:
            fb_errno = calc_forces()
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            self._forces_computed = True

    def get_charges(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        self.run_scf()
        self.charges = np.zeros(self.natoms, dtype=np.float64)
        self.shell_charges = np.zeros((self.natoms, self.nshells), dtype=np.float64)
        get_charges(self.charges, self.shell_charges.T)
        return self.charges, self.shell_charges

    def get_energies(self) -> tuple[np.float64, np.float64]:
        self.run_scf()
        self.energy, self.fermi_level = get_energies()
        return self.energy, self.fermi_level

    def get_eigenvalues(self) -> NDArray[np.float64]:
        self.run_scf()
        self.eigenvalues = np.zeros((self.nspin, self.nkpts, self.nbands))
        get_eigenvalues(self.eigenvalues[0, :, :].T)
        return self.eigenvalues

    def get_forces(self) -> NDArray[np.float64]:
        self.calc_forces()
        self.forces = np.zeros((self.natoms, 3), dtype=np.float64)
        get_forces(self.forces.T)
        return self.forces

    def update_coords(self, positions: ArrayLike) -> None:
        self.positions = np.ascontiguousarray(positions, dtype=np.float64)
        assert self.positions.shape == (self.natoms, 3), "'positions' must be a (natoms, 3) array"
        self._wrap_positions()
        if np.any(pdist(self.wrapped_positions) < 1e-6):
            raise ValueError("Atom positions are too close! (tolerance 1e-6)")
        update_coords(self.positions.T)
        self._scf_computed = False
        self._forces_computed = False
