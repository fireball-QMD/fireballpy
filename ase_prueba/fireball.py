import warnings
import numpy as np
from mendeleev import element

from ase.calculators.calculator import (Calculator,
                                        all_changes)

from _infodat import *

class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.
    """

    implemented_properties = ['energy', 'energies', 'forces']

    ignored_changes = ['initial_magmoms']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        self._numbers = set() # Keep track of our atoms

    def _generate_infodat(self, numbers):
        with open("info.dat", "w") as f:
            f.write("   dani\n")
            f.write(f"          {numbers.size:2d}  - Number of species\n")
            for i, num in enumerate(numbers):
                elem = element(num)
                f.write("======================================================================\n")
                f.write(f"  {i:2d}          - Information for this species\n")
                f.write(f"  {elem.symbol}          - Element\n")
                f.write(f"   {num:2d}          - Nuclear Z\n")
                f.write(f" {elem.atomic_weight:7.3f}          - Atomic Mass\n")
                f.write(f"  {len(shells[num]):2d}          - Number of shells; L for each shell\n")
                f.write("    " + "  ".join(shells[num]) + "\n")
                f.write(f"  {len(ps_shells[num]):2d}          - Number of shells; L for each shell  (Pseudopotential)\n")
                f.write("    " + "  ".join(ps_shells[num]) + "\n")
                f.write(f"  {cutoffs[num][0]:5.2f} - Radial cutoffs PP\n")
                f.write("   " + "   ".join(cutoffs[num][1][0]) + "\n")
                f.write("   " + "   ".join(cutoffs[num][1][1]) + "\n")
                f.write("  " + "         ".join(wffiles[num]) + "\n")
                f.write("  " + "         ".join(nafiles[num]) + "\n")
                f.write("       0.00000   - Atomic energy\n")
                f.write("======================================================================\n")

    def _calculate_energies(self, atoms):
        # Computation energy
        print("Computing energies") # Example
        self.energies = np.random.randn(len(atoms))
        self.energy = np.sum(self.energies)
        
        # Save energy
        self.results['energy'] = self.energy
        self.results['energies'] = self.energies
        self.results['free_energy'] = self.energy

    def _calculate_forces(self, atoms):
        # Requisite energies
        if 'energies' not in self.results:
            warnings.warn("Energies not computed. Computing energies", UserWarning)
            self._calculate_energies(atoms)
        
        # Computation forces
        print("Computing forces") # Example
        self.forces = np.random.randn(len(atoms))

        # Save forces
        self.results['forces'] = self.forces

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
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

    def initialize(self, atoms):
        numbers = atoms.get_atomic_numbers()
        
        # Imposes Python >=3.10
        if set(numbers) - self._numbers:
            self._numbers = set(numbers)
            self._generate_infodat(self._numbers)

        self.energies = np.empty(len(atoms))
        self.forces = np.empty((len(atoms), 3))