from __future__ import annotations
from typing import SupportsInt
import re
import numpy as np
from numpy.typing import ArrayLike

from .._errors import type_check
from ..fireball import BaseFireball

from ase import Atoms
from ase.cell import Cell
from ase.dft.kpoints import get_special_points


class BandPath:
    def __init__(self, *,
                 atoms: Atoms,
                 path: str | ArrayLike,
                 npoints: SupportsInt,
                 special_points: dict | None = None) -> None:
        # In the future we need to add new inputs and check for exclusion of them.
        # Right now we work only with ASE's Atoms.
        type_check(atoms, Atoms, 'atoms')
        fbobj = atoms.calc
        if not isinstance(fbobj, BaseFireball):
            raise RuntimeError("ASE calculator for parameter ``atoms`` must be ``Fireball``.")
        if not hasattr(fbobj, '_run_scf') or not fbobj.scf_computed:
            raise RuntimeError("SCF loop needs to be carried first")

        type_check(npoints, SupportsInt, 'npoints')
        self.npoints = int(npoints)

        # We will need an ASE cell
        self.asecell = Cell(fbobj.atomsystem.cell)
        if special_points is None:
            # Get them from ASE
            self.special_points = get_special_points(self.asecell)
        else:
            type_check(special_points, dict, 'special_points')
            self.special_points = special_points

        # See if we have a path as str or defined as an array
        if isinstance(path, str):
            # Now parse the path
            self.paths = [[name for name in re.split(r'([A-Z][a-z0-9]*)', p) if name] for p in path.split(',')]
        elif isinstance(path, np.ndarray):
            # Coordinates
            self.paths = [[f"Kpt{i}" for i in range(len(path))]]
        else:
            raise TypeError("Parameter ``path`` could not be identified.")
        kpts = self.asecell.bandpath(path, self.npoints, special_points=self.special_points).kpts  # type: ignore

        # Compute new eigenvalues
        self.fermi = fbobj.fermi_level
        self.fbobj = BaseFireball.postprocessing(fbobj, kpts=kpts)
        self.fbobj.compute_eigenvalues()
        self.eigenvalues = self.fbobj.eigenvalues
