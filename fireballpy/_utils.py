from __future__ import annotations
from ase import Atoms

from ._errors import type_check
from .fireball import BaseFireball


def fbobj_from_obj(fbobj: BaseFireball | None, atoms: Atoms | None) -> BaseFireball:
    fbobj_none = fbobj is not None
    atoms_none = atoms is not None

    nonecount = int(fbobj_none) + int(atoms_none)
    if nonecount != 1:
        raise ValueError("Only one of the Fireball-providing inputs must be specified.")

    if fbobj_none:
        type_check(fbobj, BaseFireball, 'fbobj')
        return fbobj
    if atoms_none:
        type_check(atoms, Atoms, 'atoms')
        if not hasattr(atoms, 'calc') or not isinstance(atoms.calc, BaseFireball):
            raise ValueError("ASE calculator for parameter ``atoms`` must be ``Fireball``.")
        return atoms.calc
    raise RuntimeError("Fireball-providing inputs could not be identified.")
