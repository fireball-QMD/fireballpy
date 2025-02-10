from __future__ import annotations
from copy import deepcopy
import warnings

from numpy.typing import ArrayLike
import numpy as np

from fireballpy._errors import raise_fb_error, type_check
from fireballpy._options import get_icharge, get_idipole, get_imixer
from fireballpy._correction import Correction
from fireballpy.fdata import FDataFiles
from fireballpy.atoms import AtomSystem
from fireballpy.kpoints import KPoints

from ase import Atoms

from fireballpy._fireball import (set_options,
                                  set_initial_charges,
                                  call_allocate_system,
                                  scf,
                                  get_sizes,
                                  get_initial_charges,
                                  calc_forces)

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
    atomsystem : AtomSystem
        Class with all the information of species, atomic numbers,
        positions and cell of the system.
    fdatafiles : FDataFiles
        Class with all the information on the location of the basis functions.
    kpoints : KPoints
        Class with all the information of the coordinates of the k-points as well
        as their associated weights.
    verbose : bool, optional
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
        Default is ``False``.
    total_charge : int, optional
        Total charge of the system in elementary charge units (1 extra electron would be ``total_charge = -1``).
        By default this is 0.
    correction : dict | bool | None, optional
        By default (``None``) will apply DFT3 correction for the selected
        FData if optimized parameters are available
        (see ``available_fdatas()``).
        Also, it may be a dictionary with the parameters for DFTD3 correction
        The dictionary must contain the keys ``'kind'`` (right now it can only
        be ``'dftd3'``), ``'damping'`` and either ``'method'`` or
        ``'params_tweaks'``. It may be forcefully shut down by providing ``False``.
        For more information, see
        `Simple DFT-D3 documentation <https://dftd3.readthedocs.io>`_.
    charges_method : str | None, optional
        How the autoconsistency in the charges will be performed.
        By default depends on the FData (``None``).
        If a custom FData is selected, then this parameter must be specified.
        For more information see :ref:`here <charges_methods>`.
    dipole_method : str, optional
        Whether to use (``dipole='improved'``) or not (``dipole='legacy'``)
        the improved dipole description.
        By default ``'improved'`` except for periodic systems as it is not yet implemented and thus it will be ignored.
    initial_charges : ArrayLike[float], optional
        A (natoms, nshells) array with the charges of each of the atom shells (0 if no shell for that atom is defined).
    mixer_kws : dict, optional
        Dictionary with the mixer parameters.
        More information :ref:`here <mixer>`.

    Methods
    -------
    run_scf(fix_charges=False)
        Run the SCF loop and store energy and charge related quantities.
    calc_forces(fix_charges=False)
        Compute the forces and store them.
    update_coords(positions)
        Update the coordinates in the Fortran module.

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
                 atomsystem: AtomSystem,
                 fdatafiles: FDataFiles,
                 kpoints: KPoints,
                 verbose: bool = False,
                 total_charge: int = 0,
                 correction: dict | None = None,
                 charges_method: str | None = None,
                 dipole_method: str = 'improved',
                 initial_charges: ArrayLike | None = None,
                 mixer_kws: dict | None = None) -> None:

        # Init variables
        self.scf_computed = False
        self.forces_computed = False

        # Type safety
        type_check(atomsystem, AtomSystem, 'atomsystem')
        type_check(fdatafiles, FDataFiles, 'fdatafiles')
        type_check(kpoints, KPoints, 'kpoints')
        type_check(verbose, bool, 'verbose')
        type_check(total_charge, int, 'total_charge')
        type_check(dipole_method, str, 'dipole_method')


        # Save everything to self
        self.fdatafiles = fdatafiles
        self.atomsystem = atomsystem
        self.kpoints = kpoints
        self.verbose = verbose
        self.total_charge = total_charge
        self.charges_method = charges_method
        self.dipole_method = dipole_method

        # Convenient aliases
        self.natoms = self.atomsystem.n
        self.nkpts = self.kpoints.n

        # Get charges method
        if self.charges_method is None:
            self.charges_method = self.fdatafiles.get_charges_method()
        else:
            type_check(self.charges_method, str, 'charges_method')
            self.charges_method = self.charges_method.lower()

        # Address correction
        self.correction = None
        if correction:
            type_check(correction, dict, 'correction', ' if custom parameters are desired')
            self._correction = Correction(atomsystem=self.atomsystem, **correction)
            self.correction = correction
        elif correction is None:
            self.correction = self.fdatafiles.get_correction(self.charges_method)
            if self.correction:
                self._correction = Correction(atomsystem=self.atomsystem, **self.correction)

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
        self.fdatafiles.load_fdata()
        self.atomsystem.set_coords()
        self.kpoints.set_kpoints()
        self.atomsystem.set_cell()

        # Set Fireball-like options
        self._options = {'dipole_method': np.int64(0) if self.atomsystem.isperiodic else get_idipole(dipole_method),
                         'charges_method': get_icharge(self.charges_method),
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
            self.eigenvectors = np.zeros((self.nkpts, self.norbitals, self.norbitals), dtype=np.complex128, order='C')

    def run_scf(self, fix_charges: bool = False) -> None:
        """Check if SCF is computed, execute the loop if not,
        then store energy, fermi_level, partial charges, shell charges and eigenvalues.

        Parameters
        ----------
        fix_charges : bool
            If ``True`` then only one SCF loop iteration will be performed. If ``False`` (default)
            then it will perform the usual computation.
            In general it is useful to fix the charges when making postprocessing.
        """
        if not self.scf_computed:
            converged, fb_errno, energy, fermi_level, charges = scf(self.verbose,
                                                                    fix_charges,
                                                                    self.shell_charges.T,
                                                                    self.eigenvalues.T,
                                                                    self.eigenvectors.T)
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            if not converged:
                warnings.warn("SCF loop did not converge. Try incrementing ``'max_iter'`` in ``mixer_kws``", UserWarning)
            self.energy = float(energy)
            self.fermi_level = float(fermi_level)
            self.charges = charges
            self.scf_computed = True
            if self.correction:
                self._correction.correct()  # type: ignore

    def calc_forces(self, fix_charges: bool = False) -> None:
        """Check if SCF is computed, execute the loop if not,
        then compute and store forces.

        Parameters
        ----------
        fix_charges : bool
            If ``True`` then only one SCF loop iteration will be performed. If ``False`` (default)
            then it will perform the usual computation.
            In general it is useful to fix the charges when making postprocessing.
        """
        self.run_scf(fix_charges)
        if not self.forces_computed:
            fb_errno = calc_forces(self.forces.T)
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            self.forces_computed = True

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


def fbobj_from_obj(fbobj: BaseFireball | None, atoms: Atoms | None) -> BaseFireball:
    fbobj_none = fbobj is not None
    atoms_none = atoms is not None

    nonecount = int(fbobj_none) + int(atoms_none)
    if nonecount != 1:
        raise ValueError("Only one of the Fireball-providing inputs must be specified.")

    if fbobj_none:
        type_check(fbobj, BaseFireball, 'fbobj')
        return fbobj
    if atoms_none:
        type_check(atoms, Atoms, 'atoms')
        if not hasattr(atoms, 'calc') or not isinstance(atoms.calc, BaseFireball):
            raise ValueError("ASE calculator for parameter ``atoms`` must be ``Fireball``.")
        return atoms.calc
    raise RuntimeError("Fireball-providing inputs could not be identified.")

