# TODO: bandpath return
from __future__ import annotations
from typing import Any

from numpy.typing import ArrayLike
from ase.calculators.abc import GetPropertiesMixin  # type: ignore
from ase.calculators.calculator import (Calculator, PropertyNotPresent,  # type: ignore
                                        kpts2sizeandoffsets, all_changes)
from ase.dft.kpoints import monkhorst_pack, bandpath
from ase.spectrum.band_structure import get_band_structure
import numpy as np

from .fireball import BaseFireball
from ._fdata import get_correction


class Fireball(Calculator):
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
        Name of the FData to be used. See ``available_fdatas()``
        for a table of available FData. If set to ``'custom'`` then
        an ``fdata_path`` pointing to a local FData folder must be provided.
        Please note that the first time an FData is used it needs to download
        all the necessary files.
    lazy : bool
        If set to ``True`` (default) it will only load necessary files for the species
        which are involved in the computation. If set to ``False`` it will load
        the FData for all the species available in it. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    verbose : bool
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
        Default is ``False``.
    kpts : ArrayLike[int] | dict[str, Any]
        Three element array-like structure (like list, tuple or numpy.ndarray) with three
        integer elements. A corresponding Monkhorst Pack will be generated with these indices.
        By default it will generate a (1, 1, 1) Monkhorst Pack, corresponding to the Gamma point.
        Another option is to provide a dictionary with the keys ``'path'`` and ``'npoints'`` which
        correspond to a :function:`ase.dft.kpoints.parse_path_string` parseable string and
        to the number of points in the path, respectively.
    **kwargs
        Advanced options. May requiere previous experience with DFT computations:

        +--------------------+--------------------------------------------------------------------------+
        | Property           | Description                                                              |
        +====================+==========================================================================+
        | ``gamma``          | Should the Gamma point be included?                                      |
        |                    | Yes (``True``), no (``False``), don't care (``None``, default)           |
        +--------------------+--------------------------------------------------------------------------+
        | ``charges_method`` | By default depends on the FData (``'auto'``).                            |
        |                    | For more options see :ref:`here <charges_methods>`.                      |
        +--------------------+--------------------------------------------------------------------------+
        | ``dipole_method``  | Whether to use the improved dipole description (``'improved'``, default) |
        |                    | or not (``'legacy'``). In periodic systems this is ignored as only       |
        |                    | the legacy description is implemented.                                   |
        +--------------------+--------------------------------------------------------------------------+
        | ``fdata_path``     | Path to a custom FData. Ignored unless ``fdata = 'custom'``              |
        +--------------------+--------------------------------------------------------------------------+
        | ``correction``     | By default (``'auto'``) will apply DFT3 correction for the selected      |
        |                    | FData if optimized parameters are available                              |
        |                    | (see ``available_fdatas()``). Can be turned off with ``'off'``.          |
        |                    | Also, a custom ASE calculator (usually :class:`dftd3.ase.DFTD3`) may     |
        |                    | be provided. This is preferable to using                                 |
        |                    | an :class:`ase.calculators.mixing.SumCalculator` in order to preserve    |
        |                    | the niceties of our API.                                                 |
        +--------------------+--------------------------------------------------------------------------+
        | ``mixer_kws``      | Dictionary with the charges mixer options. For reference                 |
        |                    | see :ref:`here <mixer>`.                                                 |
        +--------------------+--------------------------------------------------------------------------+

    Methods
    -------
    get_eigenvalues:
        Obtain Hamiltonian eigenvalues in eV.
    get_fermi_level:
        Obtain the Fermi energy in eV.
    get_ibz_k_points:
        Obtain the coordinates of k-points used in the computation.
    get_k_point_weights:
        Obtain the associated weights of the k-points.
    get_number_of_bands:
        Obtain the number of orbitals considered.
    get_number_of_shells:
        Obtain the number of shells considered.
    get_shell_charges:
        Obtain the charges of each shell.

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

    _prop2fun = {'energy': 'energy', 'free_energy': 'energy', 'fermi_level': 'energy',
                 'eigenvalues': 'eigenvalues', 'charges': 'charges',
                 'shell_charges': 'charges', 'forces': 'forces'}

    implemented_properties = ['energy', 'free_energy', 'forces', 'charges']

    calc_properties = ['shell_charges', 'fermi_level', 'eigenvalues',
                       'ibz_kpoints', 'kpoint_weights']

    dimensions = ['natoms', 'nspin', 'nshells', 'nspecies', 'nkpts', 'nbands']

    ignored_changes = set(['initial_magmoms'])

    def __init__(self, fdata: str, *,
                 lazy: bool = True,
                 verbose: bool = False,
                 kpts: ArrayLike | dict[str, Any] = (1, 1, 1),
                 **kwargs):

        super().__init__(**kwargs)
        assert self.parameters is not None

        # Set parameters
        self.fdata = fdata
        self.lazy = bool(lazy)
        self.verbose = bool(verbose)
        self.kpts = kpts
        self.gamma = kwargs.get('gamma', None)
        self.charges_method = kwargs.get('charges_method', 'auto')
        self.dipole_method = kwargs.get('dipole_method', 'improved')
        self.fix_charges = bool(kwargs.get('fix_charges', False))
        self.input_charges = kwargs.get('input_charges', None)
        self.mixer = kwargs.get('mixer_kws', None)
        self.fdata_path = kwargs.get('fdata_path', None)
        self._correction = kwargs.get('correction', 'auto')

        # Check kpts
        if isinstance(self.kpts, dict):
            if 'path' not in self.kpts or 'npoints' not in self.kpts:
                raise KeyError("If parameter 'kpts' is a dictionary it must have the 'path' and 'npoints' keys.")
        else:
            self.kpts = np.asarray(kpts, dtype=np.int64)
            if self.kpts.shape != (3,):
                raise ValueError("Parameter 'kpts' must be a 3-element array-like structure.")
            ksize, koffset = kpts2sizeandoffsets(self.kpts, gamma=self.gamma)
            self.kpts = monkhorst_pack(ksize) + np.array(koffset)
            self.wkpts = np.ones(self.kpts.shape[0])  # FIXED WEIGHTS

        # Check correction
        if not isinstance(self._correction, GetPropertiesMixin) and self._correction not in ['auto', 'off']:
            raise ValueError("Parameter 'correction' must be either 'auto', 'off' or an ASE Calculator")

    def _get(self, name: str) -> Any:
        if not hasattr(self, '_dft'):
            raise PropertyNotPresent(name)
        if name not in self.results:
            fname = f'_get_{self._prop2fun[name]}'
            getattr(self, fname)()
        return self.results[name]

    def _get_energy(self) -> None:
        energy, fermi_level = self._dft.get_energies()
        if self._correction is not None:
            energy += self._correct_atoms.get_potential_energy()
        self.results.update({'energy': float(energy), 'free_energy': float(energy), 'fermi_level': float(fermi_level)})

    def _get_eigenvalues(self) -> None:
        eigenvalues = self._dft.get_eigenvalues()
        self.results['eigenvalues'] = eigenvalues

    def _get_charges(self) -> None:
        charges, shell_charges = self._dft.get_charges()
        self.results.update({'charges': charges, 'shell_charges': shell_charges})

    def _get_forces(self) -> None:
        forces = self._dft.get_forces()
        if self._correction is not None:
            forces += self._correct_atoms.get_forces()
        self.results['forces'] = forces

    def _create_fireball(self) -> None:
        assert self.atoms is not None

        # Create bandpath if needed
        if isinstance(self.kpts, dict):
            self.bandpath = bandpath(cell=self.atoms.cell, **self.kpts)
            self.kpts = self.bandpath.kpts
            self.wkpts = np.ones(self.kpts.shape[0])  # FIXED WEIGHTS

        # Cell for BaseFireball interpretation
        cell = self.atoms.cell.complete().array
        for i in range(3):
            if self.atoms.pbc[i]:
                continue
            cell[i, i] = 2.0*np.abs(self.atoms.get_positions()[:, i]).max()

        self._dft = BaseFireball(
                ## FData block
                fdata=self.fdata,
                fdata_path=self.fdata_path,
                ## Atoms block
                species=set(self.atoms.get_chemical_symbols()),
                numbers=self.atoms.get_atomic_numbers(),
                positions=self.atoms.get_positions(),
                ## Periodic block
                periodic=self.atoms.pbc.any(),
                kpts=self.kpts,
                wkpts=self.wkpts,
                reducekpts=not hasattr(self, 'bandpath'),
                cell=cell,
                ## SCF block
                verbose=self.verbose,
                lazy=self.lazy,
                charges_method=self.charges_method,
                dipole_method=self.dipole_method,
                fix_charges=self.fix_charges,
                input_charges=self.input_charges,
                mixer_kws=self.mixer
        )

        # Get fixed variables
        self.nspecies, self.natoms = self._dft.nspecies, self._dft.natoms
        self.nkpts, self.nspin = self._dft.nkpts, self._dft.nspin
        self.ibz_kpoints, self.kpoint_weights = self._dft.kpts, self._dft.wkpts

        # Fetch correction if available
        if isinstance(self._correction, str) and self._correction == 'auto':
            corrtype, corrparams = get_correction(self.fdata)
            if corrtype == 'dftd3':
                try:
                    from dftd3.ase import DFTD3
                    self._correction = DFTD3(damping='d3bj', params_tweaks=corrparams)
                except ModuleNotFoundError:
                    self._correction = None
            else:
                self._correction = None
        elif isinstance(self._correction, str) and self._correction == 'off':
            self._correction = None

        # Prepare atoms for correction
        if self._correction is not None:
            self._correct_atoms = self.atoms.copy()
            # Fix DFTD3 compute volume check
            if not self.atoms.pbc.all():
                self._correct_atoms.pbc = np.array([False, False, False])
            self._correct_atoms.calc = self._correction

    def band_structure(self):
        if hasattr(self, 'bandpath'):
            return get_band_structure(calc=self, path=self.bandpath)
        return super().band_structure()

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
            If any computation has been done and eigenvalues were not yet computed.
        IndexError
            If the an out-of-bounds k-point index is requested.
        """
        spin = 0
        return self._get('eigenvalues')[spin, kpt]

    def get_fermi_level(self):
        """Get the Fermi level in electronvolts.

        Returns
        -------
        float
            Fermi energy in units of eV.

        Raises
        ------
        PropertyNotPresent
            If any computation has been done and the Fermi level was not yet computed.
        """
        return self._get('fermi_level')

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
            If any computation has been done and k-points were not yet computed.
        """
        self.results['ibz_kpoints'] = self.ibz_kpoints
        return self.ibz_kpoints

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
            If the any computation has been done and k-point weights were not yet computed.
        """
        self.results['kpoint_weights'] = self.kpoint_weights
        return self.kpoint_weights

    def get_number_of_bands(self):
        """Get the number of bands.

        Returns
        -------
        int
            Number of bands

        Raises
        ------
        PropertyNotPresent
            If the any computation has been done and the number of bands is not yet computed.
        """
        return self.get_eigenvalues().shape[2]

    def get_number_of_shells(self):
        """Get the number of shells.

        Returns
        -------
        int
            Number of shells

        Raises
        ------
        PropertyNotPresent
            If the any computation has been done and the number of shells is not yet computed.
        """
        return self.get_shell_charges().shape[1]

    def get_number_of_spins(self):
        """Get the number of spins.

        Returns
        -------
        int
            Number of spins (always 1)
        """
        return 1

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
            If the any computation has been done and the shell charges are not yet computed.
        """
        return self._get('shell_charges')

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes) -> None:
        super().calculate(atoms, properties, system_changes)
        assert atoms is not None

        # If the atoms change new BaseFireball is needed
        if 'numbers' in system_changes:
            self._create_fireball()
        # If positions change we just need to update them not repeat anything else
        elif 'positions' in system_changes:
            self._dft.update_coords(atoms.get_positions())
            if self._correction is not None:
                self._correct_atoms.positions = atoms.get_positions()

        # Compute energy always
        if 'energy' not in properties:
            properties.append('energy')
        for prop in properties:
            _ = self._get(prop)
