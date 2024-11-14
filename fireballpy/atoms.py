from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import norm, inv
from scipy.spatial.distance import pdist

from ._errors import type_check
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
    a1 : ArrayLike[float] | None
        First cell vector coordinates in angstroms.
    a2 : ArrayLike[float] | None
        Second cell vector coordinates in angstroms.
    a3 : ArrayLike[float] | None
        Third cell vector coordinates in angstroms.

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

    def __init__(self, *, species: set[str], numbers: ArrayLike, positions: ArrayLike,
                 a1: ArrayLike | None = None,
                 a2: ArrayLike | None = None,
                 a3: ArrayLike | None = None) -> None:
        # Deal with atom information
        self.sps = set(species)
        self.nums = np.ascontiguousarray(numbers, dtype=np.int64)
        if len(self.nums.shape) > 1:
            raise ValueError("Parameter ``numbers`` must be a 1D array")
        self.n = self.nums.size
        self.pos = np.ascontiguousarray(positions, dtype=np.float64)
        if self.pos.shape != (self.n, 3):
            raise ValueError("Parameter ``positions`` must be a 2D array with dimensions natoms x 3")
        # Shift coords if needed
        if self._needs_shift():
            self.pos += self._shifter

        # Deal with the cell vectors
        if a1 is None:
            a1 = [0.0, 0.0, 0.0]
        type_check(a1, ArrayLike, 'a1', ' if provided')
        if a2 is None:
            a2 = [0.0, 0.0, 0.0]
        type_check(a2, ArrayLike, 'a1', ' if provided')
        if a3 is None:
            a3 = [0.0, 0.0, 0.0]
        type_check(a3, ArrayLike, 'a1', ' if provided')
        a1 = np.ascontiguousarray(a1, dtype=np.float64)
        a2 = np.ascontiguousarray(a2, dtype=np.float64)
        a3 = np.ascontiguousarray(a3, dtype=np.float64)
        if a1.shape != (3,):
            raise ValueError("Parameter ``a1`` must have 3 elements")
        if a2.shape != (3,):
            raise ValueError("Parameter ``a1`` must have 3 elements")
        if a3.shape != (3,):
            raise ValueError("Parameter ``a1`` must have 3 elements")
        self.cell = np.ascontiguousarray([a1, a2, a3], dtype=np.float64)
        self.pbc = np.array([np.abs(c).sum() > 1e-5 for c in self.cell], dtype=np.bool_)

        # Correct cell for non-periodic dimensions
        self.isperiodic = any(self.pbc)
        for i, p in enumerate(self.pbc):
            if not p:
                self.cell[i] = np.zeros(3, dtype=np.float64, order='C')
                self.cell[i,i] = np.max(np.abs(self.pos[:,i])) + 100.0
        self.icell = 2.0*np.pi*np.ascontiguousarray(inv(self.cell).T)
        # Use the cell to check wrapped positions
        self._check_positions()

    def _needs_shift(self, tol: float = 1e-5) -> bool:
        return bool(np.any(norm(self.pos, axis=1) < tol))

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
