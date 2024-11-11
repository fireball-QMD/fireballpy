from __future__ import annotations
from itertools import product
from typing import Optional
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import norm

from _fireball import set_kpoints

class KPoints:
    """Class to hold te coordinates of the k-points and their weights.

    Parameters
    ----------
    kpts : ArrayLike[float]
        A nkpts x 3 array with the coordinates of the k-points
        in angstroms.
    wkpts : ArrayLike[float]
        A nkpts array witht the weights of each k-point.
        They do not need to be normalized.

    Methods
    -------
    reduce_kpts(tol=1e-5)
        Reduce the k-points using inverse symmetry.
        Weights are readjusted accordingly.
    set_kpoints()
        Set the k-points in the fortran module.
    """
    def __init__(self, kpts: ArrayLike, wkpts: ArrayLike) -> None:
        self.ws = np.ascontiguousarray(wkpts, dtype=np.float64)
        self.ks = np.ascontiguousarray(kpts, dtype=np.float64)
        self.n = self.ws.shape[0]
        assert self.ws.shape == (self.n,), "wkpts must be a (nkpts) size array"
        assert self.ks.shape == (self.n, 3), "kpts must be a (nkpts, 3) size array"
        self.ws /= self.ws.sum()
        self.isgamma = self.n == 1 and norm(self.ks[0, :]) < 1e-5

    @classmethod
    def from_mp(cls, mp: ArrayLike, icell: NDArray[np.float64],
                gamma: Optional[bool] = None) -> KPoints:
        """Generate KPoints from Monkhorst-Pack.

        Parameters
        ----------
        mp : ArrayLike[int]
            Three Monkhorst-Pack indices.
        icell : NDArray[float]
            A 3 x 3 array where each row corresponds to each of the reciprocal cell vectors.
        gamma : Optional[bool]
            Should the Gamma ([0, 0, 0]) point be forcefully
            included (``True``), forcefully excluded (``False``) or
            don't care whether it is included or not (``None``, default)
        """
        mpi = np.ascontiguousarray(mp, dtype=np.int64)
        assert mpi.shape == (3,), "Monkhorst-Pack must be a set of 3 positive integers"
        assert np.all(mpi > 0), "Monkhorst-Pack must be a set of 3 positive integers"
        assert icell.shape == (3, 3), "Reciprocal cell must be a (3, 3) size array"
        # Do not compute if trivial
        if mpi[0] == mpi[1] == mpi[2] == 1:
            return cls([[0.0, 0.0, 0.0]], [1.0])
        # Generate the kpoints
        kpts = np.ascontiguousarray(list(product(range(mpi[0]), range(mpi[1]), range(mpi[2]))), dtype=np.float64)
        kpts = (kpts + 0.5)/mpi
        if gamma is not None:
            kpts += 0.5/mpi * (mpi % 2 != bool(gamma))
        # Return to cartesians
        return cls(np.dot(kpts, icell), np.ones(kpts.shape[0]))

    @classmethod
    def from_density(cls, density: float, pbc: ArrayLike, icell: NDArray[np.float64],
                     gamma: Optional[bool] = None) -> KPoints:
        """Generate KPoints from a k-point density.

        Parameters
        ----------
        density : float
            K-point density in inverse angstroms.
        pbc : ArrayLike[bool]
            A boolean 3 element array specifying in which direction we
            consider periodic boundary conditions.
        icell : NDArray[float]
            A 3 x 3 array where each row corresponds to each of the reciprocal cell vectors.
        gamma : Optional[bool]
            Should the Gamma ([0, 0, 0]) point be forcefully
            included (``True``), forcefully excluded (``False``) or
            don't care whether it is included or not (``None``, default)
        """
        pbci = np.ascontiguousarray(pbc, dtype=bool)
        assert pbci.shape == (3,), "Periodic boundary conditions must be a set of 3 booleans"
        assert icell.shape == (3, 3), "Reciprocal cell must be a (3, 3) size array"
        mp = np.ones(3, dtype=np.int64)
        for i, p in enumerate(pbci):
            if not p:
                continue
            k = norm(icell[i])*density
            mp[i] = 2*int(np.ceil(0.5*k))
        return cls.from_mp(mp, icell, gamma)

    def reduce_kpts(self, tol: float = 1e-5) -> None:
        """Reduce the k-points using inverse symmetry. Weights are readjusted accordingly.

        Parameters
        ----------
        tol : float
            Tolerance for considering two k-points to be symmetric.
        """
        kpts = [self.ks[0,:]]
        wkpts = [self.ws[0]]
        for k1, w1 in zip(self.ks[1:,:], self.ws[1:]):
            for i, k2 in enumerate(kpts):
                if norm(k1 + k2) < tol:
                    wkpts[i] += w1
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


def new_kpoints(kpts: ArrayLike | float,
                wkpts: Optional[ArrayLike] = None,
                icell: Optional[NDArray[np.float64]] = None,
                pbc: Optional[ArrayLike] = None,
                gamma: Optional[bool] = None) -> KPoints:
    """Wrapper to call the appropiate KPoints constructor.

    Parameters
    ----------
    kpts : ArrayLike[int] | ArrayLike[float] | float
        Specify the k-points for a periodic computation. It can be either a set of three
        Monkhorst-Pack indices, a nkpts x 3 array with the coordinates of the k-points
        in angstroms, or a density of k-points in inverse angstroms.
    wkpts : Optional[ArrayLike[float]]
        If ``kpts`` is specified as a nkpts x 3 array with the coordinates of
        the k-points in reciprocal cell units, then ``wkpts`` may be a nkpts
        array with the weights associated to each of the k-points. By default,
        all points have the same weight.
    icell : Optional[NDArray[float]]
        A 3 x 3 array where each row corresponds to each of the reciprocal cell vectors.
        Needed if ``kpts`` is either Monkhorst-Pack indices or a density.
    pbc : Optional[ArrayLike[bool]]
        A boolean 3 element array specifying in which direction we
        consider periodic boundary conditions.
        Needed if ``kpts`` is a density.
    gamma : Optional[bool]
        Should the Gamma ([0, 0, 0]) point be forcefully
        included (``True``), forcefully excluded (``False``) or
        don't care whether it is included or not (``None``, default)

    Returns
    -------
    KPoints
        Object with the coordinates in angstroms of the k-points and their weights.

    """
    if isinstance(kpts, float):
        assert isinstance(icell, np.ndarray)
        assert pbc is not None
        return KPoints.from_density(kpts, pbc, icell, gamma)

    kpts = np.ascontiguousarray(kpts, dtype=np.float64)
    if kpts.__len__() == 3:
        assert isinstance(icell, np.ndarray)
        return KPoints.from_mp(kpts, icell, gamma)

    if wkpts is None:
        wkpts = np.ones(kpts.shape[0])
    return KPoints(kpts, wkpts)
