from __future__ import annotations
from typing import Optional
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import norm, inv
from scipy.spatial.distance import pdist

from _fireball import set_coords, set_cell, update_coords

class AtomSystem:
    """Class to contain all the information about the atomic system.

    Parameters
    ----------
    species : set[str]
        Set with the element names entering the computation.
    numbers : ArrayLike[int]
        Array with the atomic numbers of each atom.
    positions : ArrayLike[float]
        A natoms x 3 array with the positions of each atom in angstroms.
    cell : ArrayLike[float]
        A 3 x 3 array where each row is one of the lattice vectors.
    pbc : ArrayLike[bool]
        A boolean 3 element array specifying in which direction we
        consider periodic boundary conditions.

    Methods
    -------
    set_coords()
        Set the coordinate info in the Fortran module.
    set_cell()
        Set the cell info in the Fortran module.
    update_coords(positions)
        Update the coordinate info in the Fortran module.
        This is faster than setting but requires them to be setted before.
    """

    _shifter = np.ascontiguousarray([np.pi, np.exp(-1.0), np.sqrt(2.0)])

    def __init__(self, species: set[str], numbers: ArrayLike,
                 positions: ArrayLike, cell: ArrayLike, pbc: ArrayLike) -> None:
        # Save to self
        self.sps = set(species)
        self.nums = np.ascontiguousarray(numbers, dtype=np.int64)
        self.pos = np.ascontiguousarray(positions, dtype=np.float64)
        self.cell = np.ascontiguousarray(cell, dtype=np.float64)
        self.pbc = np.ascontiguousarray(pbc, dtype=np.bool_)
        self.n = self.nums.shape[0]
        # Check sizes
        assert self.nums.shape == (self.n,), "numbers must be a (nkpts) size array"
        assert self.pos.shape == (self.n, 3), "positions must be a (nkpts, 3) size array"
        assert self.cell.shape == (3, 3), "cell must be a (3, 3) size array"
        assert self.pbc.shape == (3,), "pbc must be a (3) size array"
        # Valid coords? shift if needed
        if self._needs_shift():
            self.pos += self._shifter
        # Correct cell for non-periodic dimensions
        self.isperiodic = self.pbc.any()
        for i, p in enumerate(self.pbc):
            if not p:
                self.cell[i] = np.zeros(3, dtype=np.float64, order='C')
                self.cell[i,i] = np.max(np.abs(self.pos[:,i])) + 100.0
        self.icell = 2.0*np.pi*np.ascontiguousarray(inv(self.cell).T)
        self._check_positions()

    def _needs_shift(self, tol: float = 1e-5) -> np.bool_:
        return np.any(norm(self.pos, axis=1) < tol)

    def _wrap_positions(self, eps: float = 1e-5) -> NDArray[np.float64]:
        if not self.isperiodic:
            return self.pos
        fractional = 0.5*np.dot(self.pos, self.icell)/np.pi + eps
        fractional = fractional % 1.0 - eps
        return np.dot(fractional, self.cell.T)

    def _check_positions(self, tol: float = 1e-5, eps: float = 1e-5) -> None:
        if np.any(pdist(self._wrap_positions(eps=eps)) < tol):
            raise ValueError("Atom positions are too close! (tolerance 1e-6)")

    def set_coords(self) -> None:
        """Set the coordinate info in the Fortran module.

        """
        set_coords(self.nums, self.pos.T)

    def set_cell(self) -> None:
        """Set the cell info in the Fortran module.

        """
        set_cell(self.cell[0], self.cell[1], self.cell[2])

    def update_coords(self, positions: ArrayLike) -> None:
        """Update the coordinate info in the Fortran module. This is faster than setting but requires them to be setted before.

        Parameters
        ----------
        positions : ArrayLike[float]
            A natoms x 3 array with the positions of each atom in angstroms.
        """
        self.pos = np.ascontiguousarray(positions, dtype=np.float64)
        assert self.pos.shape == (self.n, 3), "positions must be a (nkpts, 3) size array"
        self._check_positions()
        if self._needs_shift():
            self.pos += self._shifter
        update_coords(self.pos.T)


def new_atomsystem(species: set[str], numbers: ArrayLike, positions: ArrayLike,
                   a1: Optional[ArrayLike] = None,
                   a2: Optional[ArrayLike] = None,
                   a3: Optional[ArrayLike] = None) -> AtomSystem:
    """Wrapper to conveniently construct an AtomSystem object.

    Parameters
    ----------
    species : set[str]
        Set with the element names entering the computation.
    numbers : ArrayLike[int]
        Array with the atomic numbers of each atom.
    positions : ArrayLike[float]
        A natoms x 3 array with the positions of each atom in angstroms.
    a1 : Optional[ArrayLike[float]]
        First cell vector coordinates in angstroms.
    a2 : Optional[ArrayLike[float]]
        Second cell vector coordinates in angstroms.
    a3 : Optional[ArrayLike[float]]
        Third cell vector coordinates in angstroms.

    Returns
    -------
    AtomSystem
        Created AtomSystem object with all supplied information.
    """
    if a1 is None or norm(a1) < 1e-5:
        a1 = [0.0, 0.0, 0.0]
    if a2 is None or norm(a2) < 1e-5:
        a2 = [0.0, 0.0, 0.0]
    if a3 is None or norm(a3) < 1e-5:
        a3 = [0.0, 0.0, 0.0]
    a1 = np.ascontiguousarray(a1, dtype=np.float64)
    a2 = np.ascontiguousarray(a2, dtype=np.float64)
    a3 = np.ascontiguousarray(a3, dtype=np.float64)
    assert a1.shape == (3,), "a1 must have 3 elements"
    assert a2.shape == (3,), "a2 must have 3 elements"
    assert a3.shape == (3,), "a3 must have 3 elements"
    pbc = [not (a1i == a2i == a3i == 0.0) for a1i, a2i, a3i in zip(a1, a2, a3)]
    cell = np.ascontiguousarray([a1, a2, a3], dtype=np.float64)
    return AtomSystem(species, numbers, positions, cell, pbc)
