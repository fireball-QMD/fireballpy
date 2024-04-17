from __future__ import annotations
from copy import deepcopy
import errno

import numpy as np
from numpy.typing import NDArray
from ase.data import atomic_masses, chemical_symbols  # type: ignore

from ._types import integer
from fireballpy.utils import (file_as_deque,
                              read_integer_array,
                              read_line,
                              read_real_array,
                              skip_lines,
                              to_integer,
                              to_real)


class InfoDat:
    """Class to keep all the info present in the classsical info.dat file,
    needed for the adjusting the different parameters of each specie in the
    computation
    """

    def __init__(self, file_path: str):
        self.shs = {}
        self.shs_pp = {}
        self.rc_pp = {}
        self.qns = {}
        self.rcs = {}
        self.wffs = {}
        self.nafs = {}
        self.eng = {}

        dat = file_as_deque(file_path)
        skip_lines(dat)
        nsps = to_integer(read_line(dat)[0])
        for _ in range(nsps):
            skip_lines(dat, lines=3)
            z = to_integer(read_line(dat)[0])
            skip_lines(dat)
            nssh = to_integer(read_line(dat)[0])
            self.shs[z] = read_integer_array(dat, width=int(nssh))
            nssh_pp = to_integer(read_line(dat)[0])
            self.shs_pp[z] = read_integer_array(dat, width=int(nssh_pp))
            self.rc_pp[z] = to_real(read_line(dat)[0])
            self.qns[z] = read_real_array(dat, width=int(nssh))
            self.rcs[z] = read_real_array(dat, width=int(nssh))
            self.wffs[z] = read_line(dat)
            self.nafs[z] = read_line(dat)
            self.eng[z] = to_real(read_line(dat)[0])
            skip_lines(dat)

        self.anums = np.sort(list(self.shs.keys()))
        self.nsps = self.anums.size
        self.nshs = {num: self.shs[num].size for num in self.anums}
        self.nshs_pp = {num: self.shs_pp[num].size for num in self.anums}
        self.maxshs = max(
            max(self.nshs[num] for num in self.anums),
            max(self.nshs_pp[num] for num in self.anums)
        )

    def select(self, anums: NDArray[integer]) -> InfoDat:
        """Get a subset with the InfoDat with only some elements

        Parameters
        ----------
        anums: NDArray[integer]
            Array with the atomic numbers to select

        Raises:
        -------
        KeyError
            If ``anums`` contains an atomic number which is not present
        """

        uanums = np.unique(anums)
        for z in uanums:
            if z not in self.anums:
                raise KeyError(errno.EINVAL,
                               f"{z} is not present in this InfoDat")

        new_infodat = deepcopy(self)
        to_delete = set(self.anums) - set(anums)
        dicts = ("shs", "shs_pp", "rc_pp", "qns", "rcs",
                 "wffs", "nafs", "eng", "nshs", "nshs_pp")
        for attr in dicts:
            for z in to_delete:
                getattr(new_infodat, attr).pop(z)
        new_infodat.anums = uanums
        new_infodat.nsps = uanums.size
        new_infodat.maxshs = max(
            max(new_infodat.nshs[num] for num in uanums),
            max(new_infodat.nshs_pp[num] for num in uanums)
        )
        return new_infodat

    def write_ascii(self, file_path: str) -> None:
        """Write the contents of this class into a classsic
        info.dat file to be read by fireball

        Parameters
        ----------
        file_path : str
            Path to the file where the information will be saved

        Raises
        ------
        FileNotFoundError
            If ``fpath`` points to a non-valid location
        """

        info = ["  fireballpy", f"{self.nsps:12d}  - Number of species"]
        for i, z in enumerate(self.anums):
            info.extend([70*"=",
                         f"  {i+1:2d}{10*' '}- Information for this species ",
                         f"  {chemical_symbols[z]:2s}{10*' '}- Element ",
                         f"  {z:3d}{10*' '}- Nuclear Z ",
                         f"  {atomic_masses[z]:7.3f}{5*' '}- Atomic Mass ",
                         f"  {self.nshs[z]:2d}{10*' '}- Number of shells;" +
                         " L for each shell ",
                         "  " + "  ".join(map(str, self.shs[z])),
                         f"  {self.nshs_pp[z]:2d}{10*' '}- Number of shells;" +
                         " L for each shell  (Pseudopotential) ",
                         "  " + "  ".join(map(str, self.shs_pp[z])),
                         f"  {self.rc_pp[z]:5.2f} - Radial cutoffs PP ",
                         "  " + "  ".join([f"{x:5.2f}" for x in self.qns[z]]) +
                         " - Occupation numbers ",
                         "  " + "  ".join([f"{x:5.2f}" for x in self.rcs[z]]) +
                         " - Radial cutoffs ",
                         "  " + "  ".join(f"{x:<25}" for x in self.wffs[z]),
                         "  " + "  ".join(f"{x:<25}" for x in self.nafs[z]),
                         f"  {self.eng[z]:12.5f}   - Atomic energy",
                         70*"="])

        with open(file_path, "w+") as fp:
            fp.write("\n  ".join(info))
