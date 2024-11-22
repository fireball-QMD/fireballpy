from __future__ import annotations
from copy import deepcopy
from typing import SupportsFloat
import warnings

from numpy.typing import ArrayLike
import numpy as np

from ._errors import raise_fb_error, type_check
from ._options import get_icharge, get_idipole, get_imixer
from ._correction import Correction
from .fdata import FDataFiles
from .atoms import AtomSystem
from .kpoints import KPoints

from _fireball import (set_options,
                       set_initial_charges,
                       call_allocate_system,
                       scf,
                       get_sizes,
                       get_initial_charges,
                       get_charges,
                       get_energy,
                       get_eigenvalues,
                       calc_forces,
                       get_forces)

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
        Set with the element names entering the computation.
    numbers : ArrayLike[int]
        Array with the atomic numbers of each atom.
    positions : ArrayLike[float]
        A natoms x 3 array with the positions of each atom in angstroms.
    kpts : ArrayLike[int] | ArrayLike[float] | float
        Specify the k-points for a periodic computation. It can be either a set of three
        Monkhorst-Pack indices, a nkpts x 3 array with the coordinates of the k-points
        in reciprocal cell units, or a density of k-points in inverse angstroms.
        By default, it will generate a (1, 1, 1) Monkhorst-Pack corresponding to the Gamma
        point alone.
    fdata_path : str | None
        Path to a custom FData. Ignored unless ``fdata = 'custom'``.
    lazy : bool
        If set to ``True`` it will only load necessary files for the species
        which are involved in the computation. If set to ``False`` it will load
        the FData for all the species available in it. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    a1 : ArrayLike[float] | None
        First cell vector coordinates in angstroms.
    a2 : ArrayLike[float] | None
        Second cell vector coordinates in angstroms.
    a3 : ArrayLike[float] | None
        Third cell vector coordinates in angstroms.
    gamma : bool | None
        If the k-points are to be generated from a Monkhorst-Pack or from a
        k-point density, should the Gamma (``[0, 0, 0]``) point be forcefully
        included (``True``), forcefully excluded (``False``) or
        don't care whether it is included or not (``None``, default)
    verbose : bool
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
        Default is ``False``.
    total_charge : int
        Total charge of the system in elementary charge units (1 extra electron would be ``total_charge = -1``).
        By default this is 0.
    correction : dict | None
        By default (``None``) will apply DFT3 correction for the selected
        FData if optimized parameters are available
        (see ``available_fdatas()``).
        Also, it may be a dictionary with the parameters for DFTD3 correction
        The dictionary must contain the keys ``'kind'`` (right now it can only
        be ``'dftd3'``), ``'damping'`` and either ``'method'`` or
        ``'params_tweaks'``. It may be forcefully shut down by providing ``False``.
        For more information, see
        `Simple DFT-D3 documentation <https://dftd3.readthedocs.io>`_.
    charges_method : str | None
        How the autoconsistency in the charges will be performed.
        By default depends on the FData (``None``).
        If a custom FData is selected, then this parameter must be specified.
        For more information see :ref:`here <charges_methods>`.
    dipole_method : str
        Whether to use (``dipole='improved'``) or not (``dipole='legacy'``)
        the improved dipole description.
        By default ``'improved'`` except for periodic systems as it is not yet implemented and thus it will be ignored.
    fix_charges : bool
        Flag to set the initial charges as fixed for the computation, not doing SCF (default is ``False``).
    initial_charges : ArrayLike[float] | None
        A (natoms, nshells) array with the charges of each of the atom shells (0 if no shell for that atom is defined).
    mixer_kws : dict | None
        Dictionary with the mixer parameters.
        More information :ref:`here <mixer>`.

    Methods
    -------
    compute_charges()
        Get partial and shell charges.
    compute_energy()
        Get the potential energy and the fermi level in eV.
    compute_eigenvalues()
        Get the Hamiltonian eigenvalues in eV for each k-point.
    compute_forces()
        Get the forces in each atom in eV/angstroms.
    update_coords(positions)
        Update the coordinates in the Fortran module.
    postprocessing(cls, fbobj, **kwargs)
        Create a BaseFireball object with fixed charges for postprocessing.

    Notes
    -----
    If you use this code in a publication please cite us [Lewis2011]_.

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
                 kpts: ArrayLike | SupportsFloat,
                 fdata_path: str | None = None,
                 lazy: bool = True,
                 a1: ArrayLike | None = None,
                 a2: ArrayLike | None = None,
                 a3: ArrayLike | None = None,
                 gamma: bool | None = None,
                 verbose: bool = False,
                 total_charge: int = 0,
                 correction: dict | None = None,
                 charges_method: str | None = None,
                 dipole_method: str = 'improved',
                 fix_charges: bool = False,
                 initial_charges: ArrayLike | None = None,
                 mixer_kws: dict | None = None) -> None:

        # Init variables
        self.scf_computed = False
        self.forces_computed = False

        # Type safety
        type_check(lazy, bool, 'lazy')
        type_check(verbose, bool, 'verbose')
        type_check(total_charge, int, 'total_charge')
        type_check(dipole_method, str, 'dipole_method')
        type_check(fix_charges, bool, 'fix_charges')

        # Create convenient subobjects
        self.fdatafiles = FDataFiles(fdata=fdata, fdata_path=fdata_path)
        self.atomsystem = AtomSystem(species=species,
                                     numbers=numbers,
                                     positions=positions,
                                     a1=a1, a2=a2, a3=a3)
        self.kpoints = KPoints(kpts=kpts,
                               atomsystem=self.atomsystem,
                               gamma=gamma)
        self.natoms = self.atomsystem.n
        self.nkpts = self.kpoints.n

        # Save everything to self
        self.lazy = lazy
        self.verbose = verbose
        self.total_charge = total_charge
        self.charges_method = charges_method.lower()
        self.dipole_method = dipole_method
        self.fix_charges = fix_charges

        # Get charges method
        if self.charges_method is None:
            self.charges_method = self.fdatafiles.get_charges_method()
        else:
            type_check(self.charges_method, str, 'charges_method')

        # Address correction
        self.correction = None
        if correction:
            type_check(correction, dict, 'correction', ' if custom parameters are desired')
            self._correction = Correction(atomsystem=self.atomsystem, **correction)
            self.correction = correction
        elif correction is None:
            self.correction, self._correction = self.fdatafiles.get_correction(self.atomsystem, self.charges_method)

        # Address mixer
        self.mixer_kws = deepcopy(DEFAULT_MIXER)
        if mixer_kws is not None:
            type_check(mixer_kws, dict, 'mixer_kws')
            try:
                self.mixer_kws = {'method': mixer_kws['method']}
            except KeyError:
                raise ValueError("Parameter ``mixer_kws`` if specified must contain the key ``'method'``.")
            for prop in ['max_iter', 'mix_order', 'beta', 'tol', 'w0']:
                self.mixer_kws[prop] = mixer_kws.get(prop, DEFAULT_MIXER[prop])

        # Save variables to modules
        self.fdatafiles.load_fdata(self.atomsystem, self.lazy)
        self.atomsystem.set_coords()
        if not self.fix_charges:
            self.kpoints.reduce_kpts()
        self.kpoints.set_kpoints()
        self.atomsystem.set_cell()

        # Set Fireball-like options
        self._options = {'dipole_method': np.int64(0) if self.atomsystem.isperiodic else get_idipole(dipole_method),
                         'charges_method': get_icharge(self.charges_method),
                         'fix_charges': np.int64(self.fix_charges),
                         'ismolecule': np.int64(not self.atomsystem.isperiodic),
                         'isgamma': np.int64(self.kpoints.isgamma),
                         'total_charge': np.int64(-self.total_charge),
                         'mixer_method': get_imixer(self.mixer_kws['method']),
                         'max_iter': np.int64(self.mixer_kws['max_iter']),
                         'mix_order': np.int64(self.mixer_kws['mix_order']),
                         'beta': np.float64(self.mixer_kws['beta']),
                         'tol': np.float64(self.mixer_kws['tol']),
                         'w0': np.float64(self.mixer_kws['w0'])}
        set_options(**self._options)

        # Allocate module
        call_allocate_system()
        self.nbands = 0
        self._alloc_arrays()
        if initial_charges is not None:
            if not isinstance(initial_charges, (tuple, list, np.ndarray)):
                raise ValueError("Parameter ``initial_charges`` need to be coerced into an array.")
            self.initial_charges = np.ascontiguousarray(initial_charges, dtype=np.float64)
            if self.initial_charges.shape != (self.natoms, self.nshells):
                raise ValueError("Parameter ``initial_charges`` must be a natoms x nshells array. "
                                 f"Got shape {self.initial_charges.shape}")
            set_initial_charges(self.initial_charges.T)
        else:
            self.initial_charges = np.zeros((self.natoms, self.nshells), dtype=np.float64, order='C')
            get_initial_charges(self.initial_charges.T)

    def _alloc_arrays(self) -> None:
        self.nshells, self.nbands_new, self.norbitals = get_sizes()
        if self.nbands_new != self.nbands:
            self.nbands = self.nbands_new
            self.charges = np.zeros(self.natoms, dtype=np.float64, order='C')
            self.forces = np.zeros((self.natoms, 3), dtype=np.float64, order='C')
            self.shell_charges = np.zeros((self.natoms, self.nshells), dtype=np.float64, order='C')
            self.eigenvalues = np.zeros((self.nkpts, self.nbands), dtype=np.float64, order='C')

    def _run_scf(self) -> None:
        if not self.scf_computed:
            fb_errno, converged = scf(self.verbose)
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            if not converged:
                warnings.warn("SCF loop did not converge. Try incrementing ``'max_iter'`` in ``mixer_kws``", UserWarning)
            self._alloc_arrays()
            self.scf_computed = True
            if self.correction:
                self._correction.correct()  # type: ignore

    def _calc_forces(self) -> None:
        self._run_scf()
        if not self.forces_computed:
            fb_errno = calc_forces()
            if fb_errno != 0 and not self.fix_charges:
                raise_fb_error(fb_errno)
            self.forces_computed = True

    def compute_charges(self) -> None:
        """Check if SCF is computed, execute the loop if not,
        then compute partial charges and shell charges.
        """
        self._run_scf()
        get_charges(self.charges, self.shell_charges.T)

    def compute_energy(self) -> None:
        """Check if SCF is computed, execute the loop if not,
        then compute the potential energy and the fermi level.
        Automatically applies DFT-D3 correction if applicable.
        """
        self._run_scf()
        self.energy, self.fermi_level = get_energy()
        if self.correction:
            self.energy += self._correction.res['energy']  # type: ignore
        self.free_energy = self.energy

    def compute_eigenvalues(self) -> None:
        """Check if SCF is computed, execute the loop if not,
        then compute the eigenvalues.
        """
        self._run_scf()
        get_eigenvalues(self.eigenvalues.T)

    def compute_forces(self) -> None:
        """Check if SCF is computed, execute the loop if not,
        then compute forces.
        """
        self._calc_forces()
        get_forces(self.forces.T)
        if self.correction:
            self.forces -= self._correction.res['gradient']  # type: ignore

    def update_coords(self, positions: ArrayLike) -> None:
        """Update the coordinates in the module.

        Parameters
        ----------
        positions : ArrayLike[float]
            A natoms x 3 array with the new positions in angstroms.
        """
        self.atomsystem.update_coords(positions)
        if self.correction:
            self._correction.update_coords(self.atomsystem)  # type: ignore
        self.scf_computed = False
        self.forces_computed = False

    @classmethod
    def postprocessing(cls, fbobj: BaseFireball, **kwargs) -> BaseFireball:
        """Create a BaseFireball object with fixed charges for postprocessing.

        Parameters
        ----------
        fbobj : BaseFireball
            BaseFireball object to be cloned.
        **kwargs
            Parameters to substitute
        """
        type_check(fbobj, BaseFireball, 'fbobj')
        fbobj.compute_charges()

        arg_dict = {'fdata': 'custom',
                    'species': fbobj.atomsystem.species,
                    'numbers': fbobj.atomsystem.numbers,
                    'positions': fbobj.atomsystem.positions,
                    'kpts': fbobj.kpoints.kpts,
                    'fdata_path': fbobj.fdatafiles.path,
                    'lazy': fbobj.lazy,
                    'a1': fbobj.atomsystem.cell[0],
                    'a2': fbobj.atomsystem.cell[1],
                    'a3': fbobj.atomsystem.cell[2],
                    'gamma': fbobj.kpoints.gamma,
                    'verbose': fbobj.verbose,
                    'total_charge': fbobj.total_charge,
                    'correction': fbobj.correction,
                    'charges_method': fbobj.charges_method,
                    'dipole_method': fbobj.dipole_method,
                    'mixer_kws': fbobj.mixer_kws}
        for k in kwargs:
            if k not in arg_dict:
                raise ValueError(f"Parameter ``{k}`` not recognized.")
            arg_dict[k] = kwargs[k]

        return cls(fix_charges=True, initial_charges=fbobj.shell_charges, **arg_dict)
