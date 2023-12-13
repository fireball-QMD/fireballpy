import numpy as np

from ase.calculators.calculator import (Calculator,
                                        PropertyNotImplementedError,
                                        CalculatorSetupError,
                                        all_changes)

class Fireball(Calculator):
    """Python implementation of Fireball QM/MM code.
    """

    implemented_properties = ['energy', 'energies', 'forces']

    ignored_changes = ['initial_magmoms']

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)
        self._loaded_fdata = [] # We keep track of loaded Fdata's

    def _load_fdata(self):
        print(f"Loading Fdata for atoms {self._num_unique}") # Example
        self._loaded_fdata.append("fdata1.dat")

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
            
        # If the atoms change load needed Fdata and allocate memory
        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        # Check if energy needed
        if ('energy' in properties) or ('energies' in properties):
            # Computation energy
            print("Computing energies") # Example
            self.energies = np.random.randn(len(atoms))
            self.energy = np.sum(self.energies)
            
            # Save energy
            self.results['energy'] = self.energy
            self.results['energies'] = self.energies
            self.results['free_energy'] = self.energy

        # Check if forces needed
        if 'forces' in properties:
            # Requisite energies
            if 'energies' not in self.results:
                raise CalculatorSetupError("Energies must be computed before forces")
            
            # Computation forces
            print("Computing forces") # Example
            self.forces = np.random.randn(len(atoms))

            # Save forces
            self.results['forces'] = self.forces

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers()
        self._num_unique = np.unique(self.numbers)
        self._load_fdata()

        self.energies = np.empty(len(atoms))
        self.forces = np.empty((len(atoms), 3))