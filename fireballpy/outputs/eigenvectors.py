from __future__ import annotations
import numpy as np
from numpy.typing import NDArray

from fireballpy._fireball import get_orbitals
from fireballpy.fireball import BaseFireball, fbobj_from_obj
from fireballpy.outputs import OrbitalVector

from ase import Atoms


def get_eigenvectors(*, kpt: int = 0, atoms: Atoms | None = None,
                    fbobj: BaseFireball | None = None) -> NDArray:
    """Function to extract the eigenvectors from Fireball.

    This function accepts any fireballpy-like object and returns a numpy array
    with the eigenvectors for each eigenvalue as :class:`OrbitalVector` instances.
    If they are desired as a matrix, it is better to obtain them directly from the computation API.

    Parameters
    ----------
    kpt : int
        Index of the k-point to obtain the eigenvectors.
        By default takes the first (``kpt = 0``).
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
    Eigenvectors
        Numpy array with the eigenvectors as :class:`OrbitalVector`.

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
    get_orbitals(orbitals)
    eigenvectors = fbobj.eigenvectors

    return np.array([OrbitalVector(vector=ei, orbitals=orbitals) for ei in eigenvectors[kpt]], dtype=OrbitalVector)
