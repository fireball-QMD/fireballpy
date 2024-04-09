from __future__ import annotations
from dataclasses import dataclass, field

import errno

import numpy as np
from numpy.typing import NDArray
from ase.data import atomic_masses, chemical_symbols

from utils import (read_file,
                   read_float_array,
                   read_int_array,
                   read_float_entry,
                   read_int_entry,
                   read_line,
                   skip_lines)


@dataclass
class InfoDat:
    """Class to keep all the info present in the classsical info.dat file,
    needed for the adjusting the different parameters of each specie in the
    computation
    """

    shs: dict[int, np.ndarray[int, int]]
    shs_pp: dict[int, np.ndarray[int, int]]
    rc_pp: dict[int, float]
    qns: dict[int, np.ndarray[int, float]]
    rcs: dict[int, np.ndarray[int, float]]
    wffs: dict[int, list[str]]
    nafs: dict[int, list[str]]
    eng: dict[int, float]

    nsps: int = field(init=False)
    anums: list[int] = field(init=False)
    nshs: dict[int, int] = field(init=False)
    nshs_pp: dict[int, int] = field(init=False)
    maxshs: int = field(init=False)

    def __post_init__(self) -> None:
        assert self.shs.keys() \
            == self.shs_pp.keys() \
            == self.rc_pp.keys() \
            == self.qns.keys() \
            == self.rcs.keys() \
            == self.wffs.keys() \
            == self.nafs.keys() \
            == self.eng.keys(), \
            "The keys of the dictionaries are not equal"
        anums = list(self.shs.keys())
        anums.sort()
        self.anums = anums
        self.nsps = len(anums)

        nshs = {}
        nshs_pp = {}
        for num in self.anums:
            assert len(self.shs[num]) \
                == len(self.qns[num]) \
                == len(self.rcs[num]) \
                == len(self.wffs[num]) \
                == len(self.nafs[num]) - 1, \
                f"Number of items in dicts do not match for element {num}"
            nshs[num] = len(self.shs[num])
            nshs_pp[num] = len(self.shs_pp[num])
        self.nshs = nshs
        self.nshs_pp = nshs_pp

        self.maxshs = max(
            max(nshs[num] for num in self.anums),
            max(nshs_pp[num] for num in self.anums)
        )

    @classmethod
    def load(cls, fpath: str) -> InfoDat:
        new_shs = {}
        new_shs_pp = {}
        new_rc_pp = {}
        new_qns = {}
        new_rcs = {}
        new_wffs = {}
        new_nafs = {}
        new_eng = {}

        dat = read_file(fpath, header=1)
        nsps = read_int_entry(dat)
        for _ in range(nsps):
            skip_lines(dat, lines=3)
            z = read_int_entry(dat)
            skip_lines(dat)
            nssh = read_int_entry(dat)
            new_shs[z] = read_int_array(dat, width=nssh).reshape(-1,)
            nssh_pp = read_int_entry(dat)
            new_shs_pp[z] = read_int_array(dat, width=nssh_pp).reshape(-1,)
            new_rc_pp[z] = read_float_entry(dat)
            new_qns[z] = read_float_array(dat, width=nssh).reshape(-1,)
            new_rcs[z] = read_float_array(dat, width=nssh).reshape(-1,)
            new_wffs[z] = read_line(dat)
            new_nafs[z] = read_line(dat)
            new_eng[z] = read_float_entry(dat)
            skip_lines(dat)

        return cls(new_shs,
                   new_shs_pp,
                   new_rc_pp,
                   new_qns,
                   new_rcs,
                   new_wffs,
                   new_nafs,
                   new_eng)

    def select(self, anums: NDArray[int]) -> InfoDat:
        """Get a subset with the InfoDat with only some elements

        Parameters
        ----------
        anums: NDArray[int]
            Array with the atomic numbers to select

        Raises:
        -------
        KeyError
            If ``anums`` contains an atomic number which is not present
        """

        for z in anums:
            if z not in self.anums:
                raise KeyError(errno.EINVAL,
                               f"{z} is not present in this InfoDat")

        new_shs = {}
        new_shs_pp = {}
        new_rc_pp = {}
        new_qns = {}
        new_rcs = {}
        new_wffs = {}
        new_nafs = {}
        new_eng = {}
        for z in np.sort(anums):
            new_shs[z] = self.shs[z]
            new_shs_pp[z] = self.shs_pp[z]
            new_rc_pp[z] = self.rc_pp[z]
            new_qns[z] = self.qns[z]
            new_rcs[z] = self.rcs[z]
            new_wffs[z] = self.wffs[z]
            new_nafs[z] = self.nafs[z]
            new_eng[z] = self.eng[z]

        return InfoDat(new_shs,
                       new_shs_pp,
                       new_rc_pp,
                       new_qns,
                       new_rcs,
                       new_wffs,
                       new_nafs,
                       new_eng)

    def write(self, fpath: str) -> None:
        """Write the contents of this class into a classsic
        info.dat file to be read by fireball

        Parameters
        ----------
        fpath : str
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

        with open(fpath, "w+") as fp:
            fp.write("\n  ".join(info))


default_shs = {
    1: np.array([0], dtype=int),
    5: np.array([0, 1], dtype=int),
    6: np.array([0, 1], dtype=int),
    7: np.array([0, 1], dtype=int),
    8: np.array([0, 1], dtype=int),
}

default_shs_pp = {
    1: np.array([0], dtype=int),
    5: np.array([0, 1], dtype=int),
    6: np.array([0, 1], dtype=int),
    7: np.array([0, 1], dtype=int),
    8: np.array([0, 1], dtype=int),
}

default_rc_pp = {
    1: 0.2,
    5: 1.15,
    6: 0.87,
    7: 0.78,
    8: 0.71,
}

default_qns = {
    1: np.array([1.0]),
    5: np.array([2.0, 1.0]),
    6: np.array([2.0, 2.0]),
    7: np.array([2.0, 3.0]),
    8: np.array([2.0, 4.0]),
}

default_rcs = {
    1: np.array([5.42]),
    5: np.array([5.00, 5.00]),
    6: np.array([5.95, 5.95]),
    7: np.array([5.42, 5.42]),
    8: np.array([5.32, 5.32]),
}

default_wffs = {
    1: ["cinput/001_542.wf1"],
    5: ["cinput/005_500.wf1", "cinput/005_500.wf2"],
    6: ["cinput/006_595.wf1", "cinput/006_595.wf2"],
    7: ["cinput/007_542.wf1", "cinput/007_542.wf2"],
    8: ["cinput/008_532.wf1", "cinput/008_532.wf2"],
}

default_nafs = {
    1:  ["cinput/001_542.na0",
         "cinput/001_542.na1"],
    5:  ["cinput/001_500.na0",
         "cinput/001_500.na1",
         "cinput/001_500.na2"],
    6:  ["cinput/006_595.na0",
         "cinput/006_595.na1",
         "cinput/006_595.na2"],
    7:  ["cinput/007_542.na0",
         "cinput/007_542.na1",
         "cinput/007_542.na2"],
    8:  ["cinput/008_532.na0",
         "cinput/008_532.na1",
         "cinput/008_532.na2"],
}

default_eng = {
    1: 0.0,
    5: 0.0,
    6: 0.0,
    7: 0.0,
    8: 0.0,
}

default_infodat = InfoDat(default_shs,
                          default_shs_pp,
                          default_rc_pp,
                          default_qns,
                          default_rcs,
                          default_wffs,
                          default_nafs,
                          default_eng)
