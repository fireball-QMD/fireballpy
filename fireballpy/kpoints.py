from __future__ import annotations
from itertools import product
from typing import SupportsFloat
import numpy as np
from numpy.typing import ArrayLike
from numpy.linalg import norm

from ._errors import type_check
from .atoms import AtomSystem
from _fireball import set_kpoints

class KPoints:
    """Class to hold the coordinates of the k-points and their weights.

    Parameters
    ----------
    kpts : ArrayLike[int] | ArrayLike[float] | float
        Specify the k-points for a periodic computation. It can be either a set of three
        Monkhorst-Pack indices, a nkpts x 3 array with the coordinates of the k-points
        in reciprocal cell units or inverse angstroms, or a density of k-points in inverse angstroms.
    atomsystem : AtomSystem
        An AtomSystem class with the information of the atomic numbers, positions and the unit cell.
    reciprocal : bool
        If ``True`` (default), when ``kpts`` are coordinates they will be understood as reciprocal
        cell units. If ``False`` they will be understood as inverse angstroms.
        When ``kpts`` is not an array of coordinates this parameter will be ignored.
    wkpts : ArrayLike[float] | None
        If ``kpts`` is specified as a nkpts x 3 array with the coordinates of
        the k-points in reciprocal cell units, then ``wkpts`` may be a nkpts
        array with the weights associated to each of the k-points. By default,
        all points have the same weight.
    gamma : bool | None
        Should the Gamma (``[0, 0, 0]``) point be forcefully
        included (``True``), forcefully excluded (``False``) or
        don't care whether it is included or not (``None``, default)

    Methods
    -------
    reduce_kpts(tol=1e-5)
        Reduce the k-points using inverse symmetry.
        Weights are readjusted accordingly.
    set_kpoints()
        Set the k-points in the fortran module.
    """
    def __init__(self, *,
                 kpts: ArrayLike | SupportsFloat,
                 atomsystem: AtomSystem,
                 reciprocal: bool = True,
                 wkpts: ArrayLike | None = None,
                 gamma: bool | None = None) -> None:

        # Need to initialize
        mp = None
        self.isgamma = False
        self.gamma = gamma

        type_check(atomsystem, AtomSystem, 'atomsystem')
        # If kpts is an array then it's coordinates
        if isinstance(kpts, np.ndarray):
            if wkpts is None:
                wkpts = np.ones(len(kpts), dtype=np.float64)
            self.ws = np.ascontiguousarray(wkpts, dtype=np.float64)
            self.ks = np.ascontiguousarray(kpts, dtype=np.float64)
            if len(self.ws.shape) != 1:
                raise ValueError("Parameter ``wkpts`` must be a 1D array")
            self.n = self.ws.size
            if self.ks.shape != (self.n, 3):
                raise ValueError("Parameter ``kpts`` if specified as coordinates "
                                    "must be a 2D array with dimensions nkpts x 3")
            if reciprocal:
                self.rks = self.ks.copy()
                self.ks = np.dot(self.ks, atomsystem.icell)
            else:
                self.rks = 0.5*np.dot(self.ks, atomsystem.cell.T)/np.pi
            if self.n == 1 and norm(self.ks[0, :]) < 1e-5:
                self.ks = np.ascontiguousarray([[0.0, 0.0, 0.0]], dtype=np.float64)
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
        elif isinstance(kpts, SupportsFloat):
            kpts = float(kpts)
            mp = np.array([0, 0, 0], dtype=np.int64)
            for i, p in enumerate(atomsystem.pbc):
                if not p:
                    continue
                k = norm(atomsystem.icell[i])*kpts
                mp[i] = 2*np.int64(np.ceil(0.5*k))
        else:
            raise TypeError("Cannot identify the type of parameter ``kpts``.")

        # Cover the case the reading capacity of the user is non-optimal
        if mp is None:
            raise ValueError("Unrecognised type of input for ``kpts``.")

        # If it is trivial then end early
        if (mp == 1).all():
            self.ks = np.ascontiguousarray([[0.0, 0.0, 0.0]], dtype=np.float64)
            self.rks = self.ks.copy()
            self.ws = np.ascontiguousarray([1.0], dtype=np.float64)
            self.n = 1
            self.isgamma = True
            self.map = [0]
            return

        # Generate the kpoints
        kpts = np.ascontiguousarray(list(product(range(mp[0]), range(mp[1]), range(mp[2]))), dtype=np.float64)
        kpts = (kpts + 0.5)/mp
        if gamma is not None:
            kpts += 0.5/mp * (mp % 2 != bool(gamma))
        self.rks = np.ascontiguousarray(kpts, dtype=np.float64)
        self.ks = np.dot(self.rks, atomsystem.icell)
        self.ws = np.ones(self.ks.shape[0], dtype=np.float64, order='C')
        self.n = 1
        self.map = list(range(self.n))

    def reduce_kpts(self, tol: float = 1e-5) -> None:
        """Reduce the k-points using inverse symmetry. Weights are readjusted accordingly.

        Parameters
        ----------
        tol : float
            Tolerance for considering two k-points to be symmetric.
        """
        kpts = [self.ks[0]]
        wkpts = [self.ws[0]]
        for j, (k1, w1) in enumerate(zip(self.ks[1:], self.ws[1:])):
            for i, k2 in enumerate(kpts):
                if norm(k1 + k2) < tol:
                    wkpts[i] += w1
                    self.map[j + 1] = i
                    break
            else:
                kpts.append(k1)
                wkpts.append(w1)
        self.ws = np.ascontiguousarray(wkpts, dtype=np.float64)
        self.ks = np.ascontiguousarray(kpts, dtype=np.float64)
        self.n = self.ws.size

    def set_kpoints(self) -> None:
        """Set the k-points in the fortran module.
        """
        set_kpoints(self.ks.T, self.ws)
