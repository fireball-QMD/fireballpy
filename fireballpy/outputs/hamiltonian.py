from __future__ import annotations
from collections import namedtuple
import numpy as np

from fireballpy._fireball import get_orbitals, get_hs
from fireballpy.fireball import BaseFireball, fbobj_from_obj
from fireballpy.outputs import OrbitalMatrix

from ase import Atoms


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
    orbitals = np.ascontiguousarray(fbobj.natoms * [0], dtype=np.int64)
    smatt = np.zeros((fbobj.norbitals, fbobj.norbitals), dtype=np.float64, order='C')
    hmatt = np.zeros((fbobj.norbitals, fbobj.norbitals), dtype=np.float64, order='C')
    get_orbitals(orbitals)
    get_hs(smatt.T, hmatt.T)

    # Symmetrise
    hmat = 0.5*(hmatt + hmatt.T)
    smat = 0.5*(smatt + smatt.T)

    return Hamiltonian(hamiltonian=OrbitalMatrix(matrix=hmat, orbitals=orbitals),
                       overlap=OrbitalMatrix(matrix=smat, orbitals=orbitals))
