from __future__ import annotations
from typing import SupportsFloat

import numpy as np
from numpy.typing import ArrayLike
from ase.calculators.calculator import Calculator, PropertyNotPresent, all_changes

from fireballpy.atoms import AtomSystem
from fireballpy.fdata import FDataFiles
from fireballpy.kpoints import KPoints
from fireballpy.fireball import BaseFireball


class Fireball(Calculator, BaseFireball):
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
    lazy : bool, optional
        If set to ``True`` (default) it will only load necessary files for the species
        which are involved in the computation. If set to ``False`` it will load
        the FData for all the species available in it. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    verbose : bool, optional
        If ``True`` information of the convergence of the SCF loop will be printed on screen.
        Default is ``False``.
    kpts : ArrayLike[int] | ArrayLike[float] | float, optional
        Specify the k-points for a periodic computation. It can be either a set of three
        Monkhorst-Pack indices, a nkpts x 3 array with the coordinates of the k-points
        in reciprocal cell units, or a density of k-points in inverse angstroms.
        By default, it will generate a ``(1, 1, 1)`` Monkhorst-Pack corresponding to the Gamma
        point alone.
    **kwargs
        Advanced options. May requiere previous experience with DFT computations:

        +--------------------+----------------------+--------------------------------------------------------------------------+
        | Property           | Type                 | Description                                                              |
        +====================+======================+==========================================================================+
        | ``gamma``          | ``bool`` or ``None`` | If the k-points are to be generated from a Monkhorst-Pack or from a      |
        |                    |                      | k-point density, should the Gamma (``[0, 0, 0]``) point be forcefully    |
        |                    |                      | included (``True``), forcefully excluded (``False``) or                  |
        |                    |                      | don't care whether it is included or not (``None``, default)             |
        +--------------------+----------------------+--------------------------------------------------------------------------+
        | ``charges_method`` | ``str`` or ``None``  | How the autoconsistency in the charges will be performed.                |
        |                    |                      | By default depends on the FData (``None``).                              |
        |                    |                      | If a custom FData is selected, then this parameter must be specified.    |
        |                    |                      | For more options see :ref:`here <charges_methods>`.                      |
        +--------------------+----------------------+--------------------------------------------------------------------------+
        | ``dipole_method``  | ``str``              | Whether to use the improved dipole description (``'improved'``, default) |
        |                    |                      | or not (``'legacy'``). In periodic systems this is ignored as only       |
        |                    |                      | the legacy description is implemented.                                   |
        +--------------------+----------------------+--------------------------------------------------------------------------+
        | ``fdata_path``     | ``str``              | Path to a custom FData. Ignored unless ``fdata = 'custom'``.             |
        +--------------------+----------------------+--------------------------------------------------------------------------+
        | ``total_charge``   | ``int``              | Total charge of the system in elementary charge units                    |
        |                    |                      | (1 extra electron will be ``total_charge = -1``). Default is 0.          |
        +--------------------+----------------------+--------------------------------------------------------------------------+
        | ``correction``     | ``dict`` or ``None`` | By default (``None``) will apply DFT-D3 correction for the selected      |
        |                    |                      | FData if optimized parameters are available                              |
        |                    |                      | (see ``available_fdatas()``).                                            |
        |                    |                      | Also, it may be a dictionary with the parameters for DFT-D3 correction   |
        |                    |                      | The dictionary must contain the keys ``'kind'`` (right now it can only   |
        |                    |                      | be ``'dftd3'``), ``'damping'`` and either ``'method'`` or                |
        |                    |                      | ``'params_tweaks'``. It may be forcefully shut down by providing         |
        |                    |                      | ``False``. For more information, see                                     |
        |                    |                      | `Simple DFT-D3 documentation <https://dftd3.readthedocs.io>`_.           |
        +--------------------+----------------------+--------------------------------------------------------------------------+
        | ``mixer_kws``      | ``dict``             | Dictionary with the charges mixer options. For reference                 |
        |                    |                      | see :ref:`here <mixer>`.                                                 |
        +--------------------+----------------------+--------------------------------------------------------------------------+

    Methods
    -------
    get_eigenvalues(kpt=0, spin=0)
        Obtain Hamiltonian eigenvalues in eV.
    get_eigenvectors(kpt=0, spin=0)
        Obtain the wavefunction coefficients.
    get_fermi_level()
        Obtain the Fermi energy in eV.
    get_ibz_k_points()
        Obtain the coordinates of k-points used in the computation.
    get_k_point_weights()
        Obtain the associated weights of the k-points.
    get_number_of_bands()
        Obtain the number of orbitals considered.
    get_number_of_shells()
        Obtain the number of shells considered.
    get_shell_charges()
        Obtain the charges of each shell.
    fix_charges()
        Fix shell charges
    unfix_charges()
        Unfix shell charges

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
    _valid_kwargs = {'fdata_path': None, 'gamma': None, 'charges_method': None,
                     'dipole_method': 'improved', 'total_charge': 0, 'fix_charges': False,
                     'correction': None, 'initial_charges': None, 'mixer_kws': None}

    implemented_properties = ['energy', 'free_energy', 'forces', 'charges']

    calc_properties = ['shell_charges', 'fermi_level', 'eigenvalues',
                       'ibz_kpoints', 'kpoint_weights']

    ignored_changes = set(['initial_magmoms'])

    def __init__(self, fdata: str, *,
                 lazy: bool = True,
                 verbose: bool = False,
                 kpts: ArrayLike | SupportsFloat = (1, 1, 1),
                 **kwargs):
        Calculator.__init__(self, **kwargs)
        self._fb_started = False
        self.nspin = 1

        # Sadly for us we cannot check many things here, we have
        # to wait for the call to calculate.
        # We just check the most basic thing and pray the user is not mad about it.
        if fdata == 'custom' and 'fdata_path' not in kwargs:
            raise ValueError("If ``fdata='custom'`` then ``fdata_path`` must be set")
        self.fdata = fdata
        self.lazy = lazy
        self.verbose = verbose
        self.kpts = kpts
        for k in kwargs:
            if k not in self._valid_kwargs:
                raise ValueError(f"Parameter ``{k}`` not recognised.")
        self._kwargs = {}
        for k in self._valid_kwargs:
            self._kwargs[k] = kwargs.get(k, self._valid_kwargs[k])
        self.fdata_path = self._kwargs.pop('fdata_path')
        self.gamma = bool(self._kwargs.pop('gamma'))
        self.isfix_charges = bool(self._kwargs.pop('fix_charges'))

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
        if not self._fb_started:
            raise PropertyNotPresent('eigenvalues')
        return self.eigenvalues[kpt]

    def get_eigenvectors(self, kpt=0, spin=0):
        """Get the coefficients for each wave function.

        Parameters
        ----------
        kpt : int
            Index of the k-point to obtain the eigenvectors. Check ``Fireball.get_ibz_k_points()``.
            By default takes the first (``kpt = 0``).
        spin : int
            Spin index for the eigenvalues. Currently dummy as no spin polarisation is implemented.

        Returns
        -------
        NDArray[float] | NDArray[complex]
            Matrix with the eigenvectors. Each column are the coefficients of
            that wavefunction.

        Raises
        ------
        PropertyNotPresent
            If any computation has been done and eigenvectors were not yet computed.
        IndexError
            If the an out-of-bounds k-point index is requested.
        """
        if not self._fb_started:
            raise PropertyNotPresent('eigenvectors')
        return self.eigenvectors[kpt] if not self.kpoints.isgamma else np.real(self.eigenvectors[kpt])

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
        if not self._fb_started:
            raise PropertyNotPresent('fermi_level')
        return self.fermi_level

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
        if not self._fb_started:
            raise PropertyNotPresent('ibz_kpoints')
        self.results.update({'ibz_kpoints': self.kpoints.coords})
        return self.kpoints.coords

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
        if not self._fb_started:
            raise PropertyNotPresent('kpoint_weights')
        self.results.update({'kpoint_weights': self.kpoints.weights})
        return self.kpoints.weights

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
        if not self._fb_started:
            raise PropertyNotPresent('shell_charges')
        return self.nbands

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
        if not self._fb_started:
            raise PropertyNotPresent('shell_charges')
        return self.nshells

    def get_number_of_spins(self):
        """Get the number of spins.

        Returns
        -------
        int
            Number of spins (always 1)
        """
        return self.nspin

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
        if not self._fb_started:
            raise PropertyNotPresent('shell_charges')
        return self.shell_charges

    def fix_charges(self) -> None:
        """Fix shell charges for following computations
        """
        self.isfix_charges = True

    def unfix_charges(self) -> None:
        """Unfix shell charges for following computations
        """
        self.isfix_charges = False

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes) -> None:
        Calculator.calculate(self, atoms, properties, system_changes)
        assert atoms is not None

        # If the atoms change we reinit BaseFireball
        if 'numbers' in system_changes:
            atomsystem = AtomSystem(species=set(atoms.get_chemical_symbols()),
                                    numbers=atoms.get_atomic_numbers(),
                                    positions=atoms.get_positions(),
                                    a1=atoms.cell[0], a2=atoms.cell[1], a3=atoms.cell[2])
            fdatafiles = FDataFiles(fdata=self.fdata, atomsystem=atomsystem if self.lazy else None,
                                    fdata_path=self.fdata_path)
            kpoints = KPoints(kpts=self.kpts, atomsystem=atomsystem, gamma=self.gamma)
            BaseFireball.__init__(self, atomsystem=atomsystem, fdatafiles=fdatafiles, kpoints=kpoints,
                                  verbose=self.verbose, **self._kwargs)
            self._fb_started = True

        # If positions change we just need to update them not repeat anything else
        elif 'positions' in system_changes:
            self.update_coords(atoms.get_positions())

        # Compute energy, fermi_level, eigenvalues, charges and shell_charges always
        self.run_scf(self.isfix_charges)
        self.results.update({'energy': self.energy,
                             'free_energy': self.energy,
                             'fermi_level': self.fermi_level,
                             'charges': self.charges,
                             'shell_charges': self.shell_charges,
                             'eigenvalues': self.eigenvalues})

        if 'forces' in properties:
            self.calc_forces(self.isfix_charges)
            self.results.update({'forces': self.forces})
