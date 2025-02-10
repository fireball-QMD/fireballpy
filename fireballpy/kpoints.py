from __future__ import annotations
from itertools import product
from typing import SupportsFloat
import numpy as np
from numpy.typing import ArrayLike
from numpy.linalg import norm

from fireballpy._errors import type_check
from fireballpy.atoms import AtomSystem
from fireballpy._fireball import set_kpoints


class KPoints:
    """Class to hold the coordinates of the k-points and their weights.

    Parameters
    ----------
    kpts : ArrayLike[int] | ArrayLike[float] | float
        Specify the k-points for a periodic computation. It can be either a set of three
        Monkhorst-Pack indices, a nkpts x 3 array with the coordinates of the k-points
        in reciprocal cell units, or a density of k-points in inverse angstroms.
    atomsystem : AtomSystem
        An AtomSystem class with the information of the atomic numbers, positions and the unit cell.
    gamma : bool | None, optional
        Should the Gamma (``[0, 0, 0]``) point be forcefully
        included (``True``), forcefully excluded (``False``) or
        don't care whether it is included or not (``None``, default)

    Attributes
    ----------
    n : int
        Number of k-points.
    kpts : NDArray[np.float64]
        A nkpts x 3 array with the coordinates of the k-points in reciprocal cell units.
    coords : NDArray[np.float64]
        A nkpts x 3 array with the coordinates of the k-points in inverse angstroms.
    weights : NDArray[np.float64]
        Non-normalized weights associated to each k-point.
    map : list[int]
        Map from each k-point to its representative (when symmetry reduction has been performed).

    Methods
    -------
    set_kpoints()
        Set the k-points in the fortran module.
    """
    def __init__(self, *,
                 kpts: ArrayLike | SupportsFloat,
                 atomsystem: AtomSystem,
                 gamma: bool | None = None) -> None:

        # Need to initialize
        mp = None
        self.isgamma = False
        self.gamma = gamma

        type_check(atomsystem, AtomSystem, 'atomsystem')
        # If kpts is an array then it's coordinates
        if isinstance(kpts, np.ndarray):
            kpts = np.ascontiguousarray(kpts, dtype=np.float64)
            self.n = kpts.shape[0]
            if kpts.shape != (self.n, 3):
                raise ValueError("Parameter ``kpts`` if specified as coordinates "
                                 "must be a 2D array with dimensions nkpts x 3")
            self.weights = np.ones(self.n, dtype=np.float64, order='C')
            self.kpts = kpts
            self.coords = np.dot(kpts, atomsystem.icell)
            if self.n == 1 and norm(self.kpts[0, :]) < 1e-5:
                self.kpts = self.coords.copy()
                self.coords = np.ascontiguousarray([[0.0, 0.0, 0.0]], dtype=np.float64)
                self.isgamma = True
            self.map = list(range(self.n))
            return

        # If tuple or list then they are indices
        if isinstance(kpts, (tuple, list)):
            mp = np.ascontiguousarray(kpts, dtype=np.int64)
            if (mp < 1).any():
                raise ValueError("Parameter ``kpts`` if specified is Monkhorst-Pack indices "
                                 "must have values greater or equal than 1")

        # If it is a number then it should be a density
        if isinstance(kpts, SupportsFloat):
            kpts = float(kpts)
            mp = np.array([1, 1, 1], dtype=np.int64)
            for i, p in enumerate(atomsystem.pbc):
                if not p:
                    continue
                k = norm(atomsystem.icell[i])*kpts
                mp[i] = np.int64(np.ceil(k))

        # Cover the case the reading capacity of the user is non-optimal
        if mp is None:
            raise TypeError("Unrecognised type of input for ``kpts``.")

        # If it is trivial then end early
        if (mp == 1).all():
            self.coords = np.ascontiguousarray([[0.0, 0.0, 0.0]], dtype=np.float64)
            self.kpts = self.coords.copy()
            self.weights = np.ascontiguousarray([1.0], dtype=np.float64)
            self.n = 1
            self.isgamma = True
            self.map = [0]
            return

        # Generate the kpoints
        kpts = np.ascontiguousarray(list(product(range(mp[0]), range(mp[1]), range(mp[2]))), dtype=np.float64)
        kpts = (kpts + 0.5)/mp - 0.5
        if gamma is not None:
            kpts += 0.5/mp * (atomsystem.pbc * (mp % 2 != bool(gamma)))
        self.kpts = np.ascontiguousarray(kpts, dtype=np.float64)
        self.n = self.kpts.shape[0]
        self.weights = np.ones(self.n, dtype=np.float64, order='C')
        self.coords = np.dot(self.kpts, atomsystem.icell)
        self.map = list(range(self.n))
        self._reduce_kpts()

    def _reduce_kpts(self, tol: float = 1e-5) -> None:
        """Reduce the k-points using inverse symmetry. Weights are readjusted accordingly.

        Parameters
        ----------
        tol : float
            Tolerance for considering two k-points to be symmetric.
        """
        idx_unique = []
        cts_unique = []
        for i in range(self.n):
            for j in idx_unique:
                if norm(self.coords[i] + self.coords[j]) < tol:
                    self.map[i] = j
                    cts_unique[j] += 1
                    break
            else:
                idx_unique.append(i)
                cts_unique.append(1)

        self.n = len(idx_unique)
        self.weights = np.ascontiguousarray(cts_unique, dtype=np.float64)
        self.coords = np.ascontiguousarray([self.coords[i] for i in idx_unique])
        self.kpts = np.ascontiguousarray([self.kpts[i] for i in idx_unique])

    def set_kpoints(self) -> None:
        """Set the k-points in the fortran module.
        """
        set_kpoints(self.coords.T, self.weights/self.weights.sum())
