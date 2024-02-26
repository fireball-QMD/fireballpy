import warnings
from typing import Set, Optional

import numpy as np

from ase import Atoms
from ase.calculators.calculator import (Calculator,
                                        all_changes)

from .infodat import InfoDat, default_infodat

#import _fireball as _fb


class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.

    Out of the box supports molecules containing H, C, N, O and S.

    Parameters
    ----------
    infodat : InfoDat, default: default_infodat
        InfoDat object with the information of the base functions for
        the different species in the molecule. By default an InfoDat
        containing H, C, N, O and S is assumed. More information
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

    implemented_properties = ['energy', 'energies', 'forces']

    ignored_changes = ['initial_magmoms']

    def __init__(self, infodat: Optional[InfoDat] = None, **kwargs):
        Calculator.__init__(self, **kwargs)
        self._infodat = infodat if infodat is not None else default_infodat
        self._numbers = set()  # Keep track of our atoms

    # For this to work we need to remove all the info.dat stuff from
    # load_fdata and isolate the allocate in a new file called allocate.f90
    # which contains the allocation along with the deallocation if it was
    # allocated before
    def _load_infodat(self, numbers: Set[int]) -> None:
        pass
#        _fb.fdata.nspecies = len(numbers)
#        _fb.fdata.nsh_max = max(self.infodat.numshells)
#        _fb.fdata.nsh_maxPP = max(self.infodat.numshellsPP)
#        _fb.fdata.allocate()
#
#        for i, n in enumerate(numbers):
#            _fb.fdata.symbolA[i] = self.infodat.elements[n].symbol
#            _fb.fdata.nzx[i] = n
#            _fb.fdata.smass[i] = self.infodat.elements[n].atomic_weight
#            _fb.fdata.rc_PP[i] = self.infodat.cutoffPP[n]
#            _fb.fdata.etotatom[i] = self.infodat.energy[n]
#
#            _fb.fdata.nsshPP[i] = self.infodat.numshellsPP[n]
#            for j in range(self.infodat.numshellsPP[n]):
#                _fb.fdata.lsshPP[j, i] = self.infodat.shellsPP[n][j]
#            _fb.fdata.nssh[i] = self.infodat.numshells[n]
#            for j in range(self.infodat.numshells[n]):
#                _fb.fdata.lssh[j, i] = self.infodat.shells[n][j]
#                _fb.fdata.rcutoff[j, i] = self.infodat.cutoffs[n][j] * \
#                    _fb.constants.abohr
#                _fb.fdata.Qneutral[j, i] = self.infodat.qneutral[n][j]
#                _fb.fdata.wavefxn[j, i] = self.infodat.wffiles[n][j]
#                _fb.fdata.napot[j, i] = self.infodat.nafiles[n][j]
#            _fb.fdata.napot[-1, i] = self.infodat.nafiles[n][-1]
#
#        _fb.fdata.isorpmax = _fb.fdata.nssh.max()
#        _fb.fdata.isorpmax_xc = _fb.fdata.nssh.max()
#        _fb.fdata.load_fdata()

    def _calculate_energies(self, atoms: Atoms) -> None:
        # Computation energy
        print("Computing energies")  # Example
        self.energies = np.random.randn(len(atoms))
        self.energy = np.sum(self.energies)

        # Save energy
        self.results['energy'] = self.energy
        self.results['energies'] = self.energies
        self.results['free_energy'] = self.energy

    def _calculate_forces(self, atoms: Atoms) -> None:
        # Requisite energies
        if 'energies' not in self.results:
            warnings.warn("Energies not computed. Computing energies",
                          UserWarning)
            self._calculate_energies(atoms)

        # Computation forces
        print("Computing forces")  # Example
        self.forces = np.random.randn(len(atoms))

        # Save forces
        self.results['forces'] = self.forces

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes) -> None:

        Calculator.calculate(self, atoms, properties, system_changes)

        # If the atoms change load needed Fdata and allocate memory
        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        # Check if energy needed
        if ('energy' in properties) or ('energies' in properties):
            self._calculate_energies(atoms)

        # Check if forces needed
        if 'forces' in properties:
            self._calculate_forces(atoms)

    def initialize(self, atoms: Atoms) -> None:
        numbers = atoms.get_atomic_numbers()
        if set(numbers) - self._numbers:
            self._numbers = set(numbers)
            self._load_infodat(self._numbers)

        self.energies = np.empty(len(atoms))
        self.forces = np.empty((len(atoms), 3))
