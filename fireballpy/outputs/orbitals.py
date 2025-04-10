from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

from fireballpy._errors import type_check


class OrbitalVector:
    """Object to manipulate the Eigenvectors in a computation.

    This object adds the ability to select atoms and their orbitals with 2 index expressions.
    When accessing elements of the matrix one may use a ``[iatom, iorb]`` syntax.

    * ``iatom``: Refers to the index/indices of the atom.
    * ``iorb``: Refers to the orbital index/indices of the atom.

    Following :class:`np.ndarray` indexing, it is also valid to use expressions as
    ``[iatom]``, ``[iatom, iorb]``, where the absence of an index is assumed
    as a full selection.
    Slicing is supported in both the orbitals and atoms (when orbitals are not specified).

    Parameters
    ----------
    vector : NDArray
        Matrix where each row/column represents an orbital.
    orbitals : NDArray[int]
        Array with the number of orbitals each atom has.

    Attributes
    ----------
    vector : NDArray
        Matrix where each row/column represents an orbital.
    orbitals : NDArray[int]
        Array with the number of orbitals each atom has.
    """
    def __init__(self, *,
                 vector: NDArray,
                 orbitals: NDArray[np.int64]) -> None:
        type_check(vector, np.ndarray, 'NDArray')
        type_check(orbitals, np.ndarray, 'NDArray')
        self.vector = np.ascontiguousarray(vector, dtype=vector.dtype)
        if len(self.vector.shape) != 1:
            raise ValueError("Parameter ``vector`` is not a vector.")
        self.orbitals = np.ascontiguousarray(orbitals, dtype=np.int64)

        # Compute the slices for the different atoms
        self.slices = [slice(0, self.orbitals[0])]
        for i in range(1, self.orbitals.size):
            start = self.slices[-1].stop
            self.slices.append(slice(start, start + self.orbitals[i]))

    def __repr__(self) -> str:
        return self.vector.__repr__()

    def __str__(self) -> str:
        return self.vector.__str__()

    def __getitem__(self, key):
        # Only one index provided
        if not hasattr(key, '__len__'):
            key0 = []
            sls = self.slices[key]
            if isinstance(sls, slice):
                sls = [sls]
            for sl in sls:
                start = sl.start if sl.start else 0
                stop = sl.stop if sl.stop else self.vector.size
                step = sl.step if sl.step else 1
                key0.extend([i for i in range(start, stop, step)])
            return self.vector[key0]

        # Handle error
        if len(key) > 2:
            raise KeyError("Cannot interpret 3 or more indices.")

        # Atoms and orbitals
        if isinstance(key[0], slice):
            raise KeyError("Slice atom selection with orbital selection not supported. "
                           "To access concrete vector elements use the attribute ``vector``.")

        start0 = self.slices[key[0]].start
        if isinstance(key[1], slice):
            beg0 = key[1].start if key[1].start else 0
            end0 = key[1].stop if key[1].stop else self.orbitals[key[0]]
            key0 = slice(start0 + beg0, start0 + end0, key[1].step)
        else:
            key0 = start0 + key[1]

        return self.vector[key0]


class OrbitalMatrix:
    """Object to manipulate the Hamiltonian and Overlap matrix in a computation.

    This object adds the ability to select atoms and their orbitals with 4 index expressions.
    When accessing elements of the matrix one may use a ``[iatom, jatom, iorb, jorb]`` syntax.

    * ``iatom``: Refers to the index/indices of the first atom (row).
    * ``jatom``: Refers to the index/indices of the second atom (column).
    * ``iorb``: Refers to the orbital index/indices of the first atom.
    * ``jorb``: Refers to the orbital index/indices of the second atom.

    Following :class:`np.ndarray` indexing, it is also valid to use expressions as
    ``[iatom]``, ``[iatom, jatom]``, ``[iatom, jatom, iorb]``, where the absence of an index is assumed
    as a full selection.
    Slicing is supported in both the orbitals and atoms (when orbitals are not specified).

    Parameters
    ----------
    matrix : NDArray
        Matrix where each row/column represents an orbital.
    orbitals : NDArray[int]
        Array with the number of orbitals each atom has.

    Attributes
    ----------
    matrix : NDArray
        Matrix where each row/column represents an orbital.
    orbitals : NDArray[int]
        Array with the number of orbitals each atom has.
    """
    def __init__(self, *,
                 matrix: NDArray,
                 orbitals: NDArray[np.int64]) -> None:
        type_check(matrix, np.ndarray, 'NDArray')
        type_check(orbitals, np.ndarray, 'NDArray')
        self.matrix = np.ascontiguousarray(matrix, dtype=matrix.dtype)
        if len(self.matrix.shape) != 2 or (matrix.shape[0] != matrix.shape[1]):
            raise ValueError("Parameter ``matrix`` is not a matrix.")
        self.orbitals = np.ascontiguousarray(orbitals, dtype=np.int64)

        # Compute the slices for the different atoms
        self.slices = [slice(0, self.orbitals[0])]
        for i in range(1, self.orbitals.size):
            start = self.slices[-1].stop
            self.slices.append(slice(start, start + self.orbitals[i]))

    def __repr__(self) -> str:
        return self.matrix.__repr__()

    def __str__(self) -> str:
        return self.matrix.__str__()

    def __getitem__(self, key):
        # Only one index provided
        if not hasattr(key, '__len__'):
            key0 = []
            sls = self.slices[key]
            if isinstance(sls, slice):
                sls = [sls]
            for sl in sls:
                start = sl.start if sl.start else 0
                stop = sl.stop if sl.stop else self.matrix.shape[0]
                step = sl.step if sl.step else 1
                key0.extend([i for i in range(start, stop, step)])
            return self.matrix[key0, :]

        # Handle error
        if len(key) > 4:
            raise KeyError("Cannot interpret 5 or more indices.")

        # Only atoms
        if len(key) == 2:
            key0 = []
            sls = self.slices[key[0]]
            if isinstance(sls, slice):
                sls = [sls]
            for sl in sls:
                start = sl.start if sl.start else 0
                stop = sl.stop if sl.stop else self.matrix.shape[0]
                step = sl.step if sl.step else 1
                key0.extend([[i] for i in range(start, stop, step)])
            key1 = []
            sls = self.slices[key[1]]
            if isinstance(sls, slice):
                sls = [sls]
            for sl in sls:
                start = sl.start if sl.start else 0
                stop = sl.stop if sl.stop else self.matrix.shape[0]
                step = sl.step if sl.step else 1
                key1.extend([i for i in range(start, stop, step)])
            return self.matrix[key0, key1]

        # Atoms and orbitals
        if isinstance(key[0], slice) or isinstance(key[1], slice):
            raise KeyError("Slice atom selection with orbital selection not supported. "
                           "To access concrete matrix elements use the attribute ``matrix``.")

        start0 = self.slices[key[0]].start
        start1 = self.slices[key[1]].start
        if isinstance(key[2], slice):
            beg0 = key[2].start if key[2].start else 0
            end0 = key[2].stop if key[2].stop else self.orbitals[key[0]]
            step0 = key[2].step if key[2].step else 1
            key0 = slice(start0 + beg0, start0 + end0, step0)
        else:
            key0 = start0 + key[2]

        if len(key) == 4:
            if isinstance(key[3], slice):
                beg1 = key[3].start if key[3].start else 0
                end1 = key[3].stop if key[3].stop else self.orbitals[key[1]]
                step1 = key[3].step if key[3].step else 1
                key1 = slice(start1 + beg1, start1 + end1, step1)
            else:
                key1 = start1 + key[3]
        else:
            key1 = self.slices[key[1]]

        return self.matrix[key0, key1]
