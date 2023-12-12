import numpy as np

from ase.calculators.calculator import (Calculator,
                                        PropertyNotImplementedError,
                                        all_changes)

class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.
    """

    implemented_properties = ['energy', 'energies', 'forces']

    ignored_changes = ['initial_magmoms']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        self._loaded_fdata = [] # We keep track of loaded Fdata's

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers()
        self._num_unique = np.unique(self.numbers)
        self._load_fdata()

        self.energies = np.empty(len(atoms))
        self.forces = np.empty((len(atoms), 3))


    def _load_fdata(self):
        print(f"Loading Fdata for atoms {self._num_unique}")
        self._loaded_fdata.append("fdata1.dat")

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
            
        # If the atoms change or is the first time load needed Fdata and allocate memory
        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        # RIGHT NOW ASSUME ALWAYS COMPUTE ENERGY: MUST DISCUSS!!

        # Computation energy
        self.energies = np.random.randn(len(atoms)) # Example
        self.energy = np.sum(self.energies)
        
        # Save energy
        self.results['energy'] = self.energy
        self.results['energies'] = self.energies
        self.results['free_energy'] = self.energy

        # Check if forces needed
        if 'forces' in properties:
            # Computation forces
            self.forces = np.random.randn(len(atoms)) # Example

            # Save forces
            self.results['forces'] = self.forces