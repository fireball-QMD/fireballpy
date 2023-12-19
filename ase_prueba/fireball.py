import warnings
import numpy as np

from ase.calculators.calculator import (Calculator,
                                        all_changes)

class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.
    """

    implemented_properties = ['energy', 'energies', 'forces']

    ignored_changes = ['initial_magmoms']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        self._numbers = set() # Keep track of our atoms
        self._loaded_fdata = [] # We keep track of loaded Fdata's

    def _load_fdata(self, numbers):
        print(f"Loading Fdata for atoms {numbers}") # Example
        self._loaded_fdata.append("fdata1.dat")

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
        if newnums := set(numbers) - self._numbers:
            self._numbers = newnums
            self._load_fdata(newnums)

        self.energies = np.empty(len(atoms))
        self.forces = np.empty((len(atoms), 3))