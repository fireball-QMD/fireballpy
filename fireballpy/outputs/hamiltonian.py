from __future__ import annotations
from collections import namedtuple
import numpy as np
from numpy.typing import NDArray

from fireballpy._errors import type_check
from fireballpy.fireball import BaseFireball, fbobj_from_obj

from fireballpy._fireball import get_hs

from ase import Atoms


class OrbitalMatrix:
    """Object to manipulate the Hamiltonian and Overlap matrix in a computation.

    This object adds the ability to select atoms and their orbitals with 4 index expressions.
    When accessing elements of the matrix one may use a ``[iatom, jatom, iorb, jorb]`` syntax.

    * ``iatom``: Refers to the index/indices of the first atom (row).
    * ``jatom``: Refers to the index/indices of the second atom (column).
    * ``iorb``: Refers to the orbital index of the first atom.
    * ``jorb``: Refers to the orbital index of the second atom.

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
            return self.matrix[self.slices[key], :]

        # Handle error
        if len(key) > 4:
            raise KeyError("Cannot interpret 5 indices.")

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
            key0 = slice(start0 + beg0, start0 + end0, key[2].step)
        else:
            key0 = start0 + key[2]

        if len(key) == 4:
            if isinstance(key[3], slice):
                beg1 = key[3].start if key[3].start else 0
                end1 = key[3].stop if key[3].stop else self.orbitals[key[1]]
                key1 = slice(start1 + beg1, start1 + end1, key[3].step)
            else:
                key1 = start1 + key[3]
        else:
            key1 = self.slices[key[1]]

        return self.matrix[key0, key1]


Hamiltonian = namedtuple('Hamiltonian', ['hamiltonian', 'overlap'])


def get_hamiltonian(*, atoms: Atoms | None = None,
                    fbobj: BaseFireball | None = None) -> Hamiltonian:
    """Function to extract the hamiltonian and overlap matrices from Fireball.

    This function accepts any fireballpy-like object and returns a named tuple
    with the hamiltonian and the overlap matrices as :class:`OrbitalMatrix` instances.

    Parameters
    ----------
    atoms : ase.Atoms | None
        An :class:`ase.Atoms` object with a :class:`Fireball` attached as calculator.
        The SCF loop must be computed before.
        This parameter is not compatible with ``fbobj``.
    fbobj : BaseFireball | None
        Directly a low-level :class:`BaseFireball` object.
        In general this is reserved only for advanced users.
        The SCF loop must be computed before.
        This parameter is not compatible with ``atoms``.

    Returns
    -------
    Hamiltonian
        Named tuple with the hamiltonian (``.hamiltonian``) and overlap (``.overlap``) matrices.

    Raises
    ------
    ValueError
        If more than one fireballpy object is passed.
    RuntimeError
        If the SCF loop was not computed.
    """

    fbobj = fbobj_from_obj(fbobj=fbobj, atoms=atoms)
    if not hasattr(fbobj, 'run_scf') or not fbobj.scf_computed:
        raise RuntimeError("SCF loop needs to be carried first")
    fbobj.run_scf(fix_charges=True)

    # Get from module
    orbitals = np.zeros(fbobj.natoms, dtype=np.int64)
    smat = np.zeros((fbobj.norbitals, fbobj.norbitals), dtype=np.float64, order='C')
    hmat = np.zeros((fbobj.norbitals, fbobj.norbitals), dtype=np.float64, order='C')
    get_hs(orbitals, smat.T, hmat.T)

    # Symmetrise
    hmat = 0.5*(hmat + hmat.T)
    smat = 0.5*(smat + smat.T)

    return Hamiltonian(hamiltonian=OrbitalMatrix(matrix=hmat, orbitals=orbitals),
                       overlap=OrbitalMatrix(matrix=smat, orbitals=orbitals))
