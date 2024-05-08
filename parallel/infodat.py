from __future__ import annotations
from dataclasses import dataclass, field
import os

import numpy as np
from numpy.typing import NDArray
from ase.data import atomic_masses, chemical_symbols  # type: ignore
from ase.units import Bohr  # type: ignore

from ._types import integer, real
from fireballpy.utils import (alloc_integer,
                              alloc_real,
                              file_as_deque,
                              read_integer_array,
                              read_line,
                              read_real_array,
                              skip_lines,
                              to_integer,
                              to_real)

from _fireball import calc_me_max


@dataclass
class InfoDat:
    """Class to keep all the info present in the classical info.dat file,
    needed for the adjusting the different parameters of each specie in the
    computation
    """

    fdata_folder: str
    nzx: NDArray[integer]

    # info.dat related variables
    nsh_max: integer = field(init=False)
    nspecies: integer = field(init=False)
    nssh: NDArray[integer] = field(init=False)
    nsshpp: NDArray[integer] = field(init=False)
    lssh: NDArray[integer] = field(init=False)
    lsshpp: NDArray[integer] = field(init=False)
    etotatom: NDArray[real] = field(init=False)
    smass: NDArray[real] = field(init=False)
    rc_pp: NDArray[real] = field(init=False)
    rcutoff: NDArray[real] = field(init=False)
    qneutral: NDArray[real] = field(init=False)
    symbola: NDArray[np.str_] = field(init=False)
    wavefxn: NDArray[np.str_] = field(init=False)
    napot: NDArray[np.str_] = field(init=False)

    # me_max variables
    me2c_max: integer = field(init=False)
    me2cdipx_max: integer = field(init=False)
    me2cdipy_max: integer = field(init=False)
    me2cpp_max: integer = field(init=False)
    mes_max: integer = field(init=False)
    me3c_max: integer = field(init=False)

    def __post_init__(self):
        self._ensure_path()
        self._get_scalars()
        self._init_arrays()
        self._get_arrays()
        self._get_me_max()

    def _ensure_path(self):
        self.fdata_folder = os.path.expanduser(self.fdata_folder)
        if self.fdata_folder[-1] != os.path.sep:
            self.fdata_folder += os.pathsep

    def _get_scalars(self):
        self.nspecies = self.nzx.size

        dat = file_as_deque(os.path.join(self.fdata_folder, 'info.dat'))
        skip_lines(dat)
        nspecies = to_integer(read_line(dat)[0])

        nsh_max = 0
        for _ in range(nspecies):
            skip_lines(dat, 3)
            z = to_integer(read_line(dat)[0])
            if z not in self.nzx:
                skip_lines(dat, 12)
                continue
            skip_lines(dat)
            nsh_max = max(nsh_max, to_integer(read_line(dat)[0]))
            skip_lines(dat)
            nsh_max = max(nsh_max, to_integer(read_line(dat)[0]))
            skip_lines(dat, 8)
        self.nsh_max = nsh_max

    def _init_arrays(self):
        self.nssh = alloc_integer(self.nsh_max)
        self.nsshpp = alloc_integer(self.nsh_max)
        self.lssh = alloc_integer(self.nsh_max, self.nspecies)
        self.lsshpp = alloc_integer(self.nsh_max, self.nspecies)
        self.etotatom = alloc_real(self.nspecies)
        self.smass = alloc_real(self.nspecies)
        self.rc_pp = alloc_real(self.nspecies)
        self.rcutoff = alloc_real(self.nspecies, self.nsh_max)
        self.qneutral = alloc_real(self.nsh_max, self.nspecies)
        self.symbola = np.empty(self.nspecies, np.dtype('<U2'))
        self.wavefxn = np.empty((self.nsh_max, self.nspecies),
                                np.dtype('<U25'))
        self.napot = np.empty((self.nsh_max + 1, self.nspecies),
                              np.dtype('<U25'))

    def _get_arrays(self):
        dat = file_as_deque(os.path.join(self.fdata_folder, 'info.dat'))
        skip_lines(dat)
        nspecies = to_integer(read_line(dat)[0])

        i = 0
        for _ in range(nspecies):
            skip_lines(dat, 3)
            z = to_integer(read_line(dat)[0])
            if z not in self.nzx:
                skip_lines(dat, 12)
                continue
            self.symbola[i] = chemical_symbols[z].ljust(2)
            self.smass[i] = atomic_masses[z]
            skip_lines(dat)
            w = to_integer(read_line(dat)[0])
            self.nssh[i] = w
            self.lssh[:int(w), i] = read_integer_array(dat, width=int(w))
            wpp = to_integer(read_line(dat)[0])
            self.nsshpp[i] = wpp
            self.lsshpp[:int(wpp), i] = read_integer_array(dat, width=int(wpp))
            self.rc_pp[i] = to_real(read_line(dat)[0])
            self.qneutral[:int(w), i] = read_real_array(dat, width=int(w))
            self.rcutoff[i, :int(w)] = read_real_array(dat, width=int(w))*Bohr
            wavefxn = read_line(dat)
            for j in range(int(w)):
                self.wavefxn[j, i] = wavefxn[j].strip().ljust(25)
            napot = read_line(dat)
            for j in range(int(w) + 1):
                self.napot[j, i] = napot[j].strip().ljust(25)
            self.etotatom[i] = to_real(read_line(dat)[0])
            skip_lines(dat)

    def _get_me_max(self):
        res = calc_me_max(self.nssh, self.nsshpp, self.lssh, self.lsshpps)
        self.me2c_max = res[0]
        self.me2cdipx_max = res[1]
        self.me2cdipx_max = res[2]
        self.me2cpp_max = res[3]
        self.me3c_max = res[4]
        self.mes_max = res[5]
