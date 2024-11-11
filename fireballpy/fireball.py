from __future__ import annotations
from copy import deepcopy
from typing import Optional, Any

from numpy.typing import ArrayLike, NDArray
import numpy as np

from ._errors import raise_fb_error
from ._options import get_icharge, get_idipole, get_imixer
from ._correction import new_correction
from .fdata import FDataFiles
from .atoms import AtomSystem
from .kpoints import KPoints

from _fireball import (set_options,
                       set_charges,
                       call_allocate_system,
                       scf,
                       get_sizes,
                       get_charges,
                       get_energies,
                       get_eigenvalues,
                       calc_forces,
                       get_forces)

# Could be replaced by constants in future versions if
# ase dependency is to be removed
from ase.units import Hartree, Bohr

DEFAULT_MIXER = {'mixer_method': 'johnson',
                 'max_iter': np.int64(200),
                 'mix_order': np.int64(6),
                 'beta': np.float64(0.1),
                 'tol': np.float64(1e-8),
                 'w0': np.float64(0.01)}


def _type_check(obj: Any, types: tuple, name: str) -> None:
    if not isinstance(obj, types):
        raise ValueError(f"Parameter {name} must be of types {','.join(map(str, types))}")


class BaseFireball:
    """Base Fireball class to be interfaced with any high-level package such as ASE.

    Parameters
    ----------
    fdata : FDataFiles
        A FDataFiles object with the information of the location,
        species and (if provided by us) optimal parameters.
    atomsystem : AtomSystem
        An AtomSystem object with the information of the species,
        the atomic numbers, the positions and the unit cell.
    kpts: KPoints
        A KPoints object with the information of the k-point coordinates
        in reciprocal cell units and their associated weights.
    lazy : Optional[bool]
        If set to ``True`` it will only load necessary files for the species
        which are involved in the computation. If set to ``False`` it will load
        the FData for all the species available in it. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    verbose : Optional[bool]
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
        Default is ``False``.
    total_charge : Optional[int]
        Total charge of the system in elementary charge units (1 extra electron will be ``total_charge = -1``).
        By default this is 0.
    correction : Optional[str | dict]
        By default (``'auto'``) will apply DFT3 correction for the selected
        FData if optimized parameters are available
        (see ``available_fdatas()``). Can be turned off with ``'off'``.
        Also, it may be a dictionary with the parameters for DFTD3 correction
        The dictionary must contain the keys ``'type'`` (right now it can only
        be ``'dftd3'``), ``'damping'`` and either ``'method'`` or
        ``'params_tweaks'``. For more information, see
        `Simple DFT-D3 documentation <dftd3.readthedocs.io>`_.
    charges_method : Optional[str]
        How the autoconsistency in the charges will be performed.
        By default depends on the FData (``'auto'``).
        If a custom FData is selected, then this parameter must be specified.
        For more information see :ref:`here <charges_methods>`.
    dipole_method : Optional[str]
        Whether to use (``dipole='improved'``) or not (``dipole='legacy'``)
        the improved dipole description.
        By default ``'improved'`` except for periodic systems as it is not yet implemented and thus it will be ignored.
    fix_charges : Optional[bool]
        Flag to set the initial charges as fixed for the computation, not doing SCF (default is ``False``).
    initial_charges : Optional[NDArray[float]]
        A (natoms, nshells) array with the charges of each of the atom shells (0 if no shell for that atom is defined).
    mixer_kws : Optional[dict]
        Dictionary with the mixer parameters.
        More information :ref:`here <mixer>`.

    Methods
    -------
    compute_charges()
        Get partial and shell charges.
    compute_energies()
        Get the potential energy and the fermi level in eV.
    compute_eigenvalues()
        Get the Hamiltonian eigenvalues in eV for each k-point.
    compute_forces()
        Get the forces in each atom in eV/angstroms.
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
                 # Mandatory
                 fdata: FDataFiles,
                 atomsystem: AtomSystem,
                 kpts: KPoints,
                 # Optional (has default)
                 lazy: bool = True,
                 verbose: bool = False,
                 total_charge: int = 0,
                 correction: str | dict[str, Any] = 'auto',
                 charges_method: str = 'auto',
                 dipole_method: str = 'improved',
                 fix_charges: bool = False,
                 # Optional (no default)
                 initial_charges: Optional[NDArray[np.float64]] = None,
                 mixer_kws: Optional[dict[str, Any]] = None) -> None:

        # Init variables
        self._scf_computed = False
        self._forces_computed = False
        self._has_correction = False
        self.nspin = 1
        self.mixer = deepcopy(DEFAULT_MIXER)

        # Type safety
        _type_check(fdata, (FDataFiles,), 'fdata')
        _type_check(atomsystem, (AtomSystem,), 'atomsystem')
        _type_check(kpts, (KPoints,), 'kpts')
        _type_check(lazy, (bool, np.bool_), 'lazy')
        _type_check(verbose, (bool, np.bool_), 'verbose')
        _type_check(total_charge, (int, np.int_), 'total_charge')
        _type_check(correction, (str, dict), 'correction')
        _type_check(charges_method, (str,), 'charges_method')
        _type_check(dipole_method, (str,), 'dipole_method')
        _type_check(fix_charges, (bool, np.bool_), 'fix_charges')

        # Save to self
        self.fdata = fdata
        self.atomsystem = atomsystem
        self.kpts = kpts
        self.islazy = lazy
        self.verbose = verbose
        self.total_charge = total_charge
        self.correction = correction
        self.charges_method = charges_method
        self.dipole_method = dipole_method
        self.fix_charges = fix_charges

        # Address mixer
        if mixer_kws is not None:
            mixer = deepcopy(mixer_kws)
            assert isinstance(mixer, dict), "'mixer_kws' must be a dictionary"
            assert 'method' in mixer, "'mixer' must have a 'method' key"
            self.mixer = {'mixer_method': mixer['method']}
            for prop in ['max_iter', 'mix_order']:
                self.mixer[prop] = np.int64(mixer.get(prop, DEFAULT_MIXER[prop]))
            for prop in ['beta', 'tol', 'w0']:
                self.mixer[prop] = np.float64(mixer.get(prop, DEFAULT_MIXER[prop]))

        # Address correction
        if isinstance(self.correction, str):
            if self.correction not in ['auto', 'off']:
                raise ValueError("'correction' must be either 'auto', "
                                "'off' or a parameter dictionary.")
            if self.correction == 'auto':
                self.corr = self.fdata.get_correction(self.atomsystem)
                self._has_correction = self.corr is not None
        elif isinstance(self.correction, dict):
            for p in ['type', 'damping']:
                if p not in self.correction:
                    raise ValueError(f"{p} must be present in 'correction'")
            self.corr = new_correction(self.atomsystem,
                                       self.correction['type'],
                                       self.correction['damping'],
                                       method=self.correction.get('method', None),
                                       params_tweaks=self.correction.get('params_tweaks', None))
            self._has_correction = True

        # Set Fireball-like options
        if self.charges_method == 'auto':
            self.charges_method = self.fdata.get_charges_method()

        self._options = {'dipole_method': np.int64(0) if self.atomsystem.isperiodic else get_idipole(dipole_method),
                         'charges_method': get_icharge(self.charges_method),
                         'fix_charges': np.int64(self.fix_charges),
                         'ismolecule': np.int64(not self.atomsystem.isperiodic),
                         'isgamma': np.int64(self.kpts.isgamma),
                         'total_charge': np.int64(-self.total_charge)}
        self._options.update(self.mixer)
        self._options['mixer_method'] = get_imixer(self.mixer['mixer_method'])
        set_options(**self._options)

        # Save variables to modules
        self.fdata.load_fdata(self.atomsystem.sps, self.islazy)
        self.atomsystem.set_coords()
        self.atomsystem.set_cell()
        if not self.fix_charges:
            self.kpts.reduce_kpts()
        self.kpts.set_kpoints()

        # Allocate module
        call_allocate_system()
        self.nbands = 0
        self._alloc_arrays()
        if initial_charges is not None:
            self.initial_charges = np.ascontiguousarray(initial_charges, dtype=np.float64)
            assert self.initial_charges.shape == (self.atomsystem.n, self.nshells), "'initial_charges must be a (natoms, nshells) array'"
            set_charges(self.initial_charges.T)

    def _alloc_arrays(self) -> None:
        self.nshells, self.nbands_new = get_sizes()
        if self.nbands_new != self.nbands:
            self.nbands = self.nbands_new
            self.charges = np.zeros(self.atomsystem.n, dtype=np.float64, order='C')
            self.forces = np.zeros((self.atomsystem.n, 3), dtype=np.float64, order='C')
            self.shell_charges = np.zeros((self.atomsystem.n, self.nshells), dtype=np.float64, order='C')
            self.eigenvalues = np.zeros((self.nspin, self.kpts.n, self.nbands), dtype=np.float64, order='C')

    def _run_scf(self) -> None:
        if not self._scf_computed:
            fb_errno = scf(self.verbose)
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            self._alloc_arrays()
            self._scf_computed = True
            if self._has_correction:
                self.corrres = self.corr.disp.get_dispersion(param=self.corr.dpar, grad=True)  # type: ignore

    def _calc_forces(self) -> None:
        self._run_scf()
        if not self._forces_computed:
            fb_errno = calc_forces()
            if fb_errno != 0:
                raise_fb_error(fb_errno)
            self._forces_computed = True

    def compute_charges(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Check if SCF is computed, execute the loop if not,
        then return the partial charges and the shell charges.

        Returns
        -------
        charges : NDArray[np.float64]
            A natoms array with the partial charges.
        shell_charges : NDArray[np.float64]
            A natoms x nshells array with the charges for each shell.
            If a shell is not defined for some atom, its charge is set to 0.
        """
        self._run_scf()
        get_charges(self.charges, self.shell_charges.T)
        return self.charges, self.shell_charges

    def compute_energies(self) -> tuple[np.float64, np.float64]:
        """Check if SCF is computed, execute the loop if not,
        then return the potential energy and the fermi level.
        Automatically applies DFT-D3 correction if applicable.

        Returns
        -------
        energy : np.float64
            Potential energy in eV.
        fermi_level : np.float64
            Fermi level in eV.
        """
        self._run_scf()
        self.energy, self.fermi_level = get_energies()
        if self._has_correction:
            self.energy += self.corrres['energy']*Hartree
        return self.energy, self.fermi_level

    def compute_eigenvalues(self) -> NDArray[np.float64]:
        """Check if SCF is computed, execute the loop if not,
        then return the eigenvalues of the Hamiltonian.

        Returns
        -------
        eigenvalues : NDArray[np.float64]
            A 1 x nkpts x nbands array with the eigenvalues of the Hamiltonian in eV.
        """
        self._run_scf()
        get_eigenvalues(self.eigenvalues[0, :, :].T)
        return self.eigenvalues

    def compute_forces(self) -> NDArray[np.float64]:
        """Check if SCF is computed, execute the loop if not,
        then compute forces if they have not been computed before.

        Returns
        -------
        forces : NDArray[np.float64]
            A natoms x 3 array with the forces in eV/angstroms.
        """
        self._calc_forces()
        get_forces(self.forces.T)
        if self._has_correction:
            self.forces -= self.corrres['gradient']*Hartree/Bohr
        return self.forces

    def update_coords(self, positions: ArrayLike) -> None:
        """Update the coordinates in the module.

        Parameters
        ----------
        positions : ArrayLike[float]
            A natoms x 3 array with the new positions in angstroms.
        """
        self.atomsystem.update_coords(positions)
        if self._has_correction:
            self.corr.disp.update(self.atomsystem.pos/Bohr, self.atomsystem.cell/Bohr)  # type: ignore
