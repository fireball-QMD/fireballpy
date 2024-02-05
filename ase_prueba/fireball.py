import warnings
from typing import Set

import numpy as np
from mendeleev import element

from ase import Atoms
from ase.calculators.calculator import (Calculator,
                                        all_changes)

from _infodat import *

import _fireball


class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.
    """

    implemented_properties = ['energy', 'energies', 'forces']

    ignored_changes = ['initial_magmoms']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        self._numbers = set()  # Keep track of our atoms

    # For this to work we need to remove all the info.dat stuff from
    # load_fdata and isolate the allocate in a new file called allocate.f90
    # which contains the allocation along with the deallocation if it was
    # allocated before
    def _load_infodat(self, numbers: Set[int]) -> None:
        _fireball.fdata.nspecies = len(numbers)
        _fireball.fdata.nsh_max = max(len(_shells[n]) for n in numbers)
        _fireball.fdata.nsh_maxPP = max(len(_ps_shells[n]) for n in numbers)
        _fireball.fdata.allocate()

        for i, n in enumerate(numbers):
            el = element(n)
            _fireball.fdata.symbolA[i] = el.symbol
            _fireball.fdata.nzx[i] = el.atomic_number
            _fireball.fdata.smass[i] = el.atomic_weight
            _fireball.fdata.rc_PP[i] = _cutoffs[n][0]
            _fireball.fdata.etotatom[i] = 0.0

            _fireball.fdata.nsshPP[i] = len(_ps_shells[n])
            for j in range(_fireball.fdata.nsshPP[i]):
                _fireball.fdata.lsshPP[j, i] = _ps_shells[n][j]
            _fireball.fdata.nssh[i] = len(_shells[n])
            for j in range(_fireball.fdata.nssh[i]):
                _fireball.fdata.lssh[j, i] = _shells[n][j]
                _fireball.fdata.Qneutral[j, i] = _cutoffs[n][1][0, j]
                _fireball.fdata.rcutoff_temp[j, i] = _cutoffs[n][1][1, j] * \
                    _fireball.constants.abohr
                _fireball.fdata.wavefxn[j, i] = _wffiles[n][j]
                _fireball.fdata.napot[j, i] = _nafiles[n][j]
            _fireball.fdata.napot[-1, i] = _nafiles[n][-1]

        _fireball.fdata.isorpmax = _fireball.fdata.nssh.max()
        _fireball.fdata.isorpmax_xc = _fireball.fdata.nssh.max()
        _fireball.fdata.load_fdata()

#    def _generate_infodat(self, numbers):
#        with open("info.dat", "w") as f:
#            f.write("   dani\n")
#            f.write(f"          {numbers.size:2d}  - Number of species\n")
#            for i, num in enumerate(numbers):
#                elem = element(num)
#                f.write("======================================================================\n")
#                f.write(f"  {i:2d}          - Information for this species\n")
#                f.write(f"  {elem.symbol}          - Element\n")
#                f.write(f"   {num:2d}          - Nuclear Z\n")
#                f.write(f" {elem.atomic_weight:7.3f}          - Atomic Mass\n")
#                f.write(f"  {len(shells[num]):2d}          - Number of shells; L for each shell\n")
#                f.write("    " + "  ".join(shells[num]) + "\n")
#                f.write(f"  {len(ps_shells[num]):2d}          - Number of shells; L for each shell  (Pseudopotential)\n")
#                f.write("    " + "  ".join(ps_shells[num]) + "\n")
#                f.write(f"  {cutoffs[num][0]:5.2f} - Radial cutoffs PP\n")
#                f.write("   " + "   ".join(cutoffs[num][1][0]) + "\n")
#                f.write("   " + "   ".join(cutoffs[num][1][1]) + "\n")
#                f.write("  " + "         ".join(wffiles[num]) + "\n")
#                f.write("  " + "         ".join(nafiles[num]) + "\n")
#                f.write("       0.00000   - Atomic energy\n")
#                f.write("======================================================================\n")

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
