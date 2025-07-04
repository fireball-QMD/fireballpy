from __future__ import annotations
from copy import deepcopy
import os
import sys
from typing import Any, Iterable

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import numpy as np
from scipy.interpolate import CubicSpline
from ase.units import Bohr

from fireballpy import __version__ as __fb_version__
from fireballpy._errors import type_check
from fireballpy.utils import get_data_from_url, read_wf, read_wf_info, ANGULAR_MOMENTUM_REV
from fireballpy.atoms import AtomSystem
from fireballpy._fireball import loadfdata_from_path

_loaded_fdata: tuple[str | None, set[str]] = ('', set(''))

#######################################
# UTIL FUNCTIONS FOR DOWNLOADING FDATAS
def _load_fdatas() -> dict[str, Any]:
    fdata_db = os.path.join(os.path.dirname(__file__), 'fdata.toml')
    with open(fdata_db, 'rb') as fp:
        return tomllib.load(fp)
#######################################

######################################
# UTIL FUNCTIONS FOR CHECKING INFO.DAT
def _get_blocks_infodat(infodat: list[str]) -> Iterable[tuple[int, int]]:
    eqsplits = [i for i, l in enumerate(infodat) if l == 70*'=']
    return list(zip(eqsplits[:-1:2], eqsplits[1::2]))


def _get_specie_block(infodat: list[str], block: tuple[int, int]) -> str:
    return infodat[block[0]+2].split('-')[0].strip()
######################################

def available_fdatas() -> list[str]:
    fdatas = _load_fdatas()
    keys = list(fdatas.keys())
    for i, k in enumerate(keys):
        if 'correction' in fdatas[k]:
            keys[i] += f' ({fdatas[k]["correction"]["kind"].upper()} available)'
    return keys

def get_fdata(name: str) -> dict[str, Any]:
    fdatas = _load_fdatas()
    try:
        return fdatas[name]
    except KeyError:
        raise ValueError("Provided FData 'name' is not in the database")


class FDataFiles:
    """Class to contain needed information about FData files in disk

    Parameters
    ----------
    fdata : str
        Name of the FData to be used. See ``available_fdatas()``
        for a table of available FData. If set to ``'custom'`` then
        an ``fdata_path`` pointing to a local FData folder must be provided.
        Please note that the first time an FData is used it needs to download
        all the necessary files.
    atomsystem : AtomSystem | None, optional
        AtomSystem object with the information of the species which intervine in the computation.
        If ``None`` (default) then the whole FData will be loaded. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.
    fdata_path : str | None, optional
        Path to a custom FData. Ignored unless ``fdata = 'custom'``.

    Attributes
    ----------
    name : str
        Name of the provided FData.
    path : str
        Path to the FData.

    Methods
    -------
    load_fdata()
        Set Fortran module variables and load the FData into memory.
    get_charges_method()
        Return the string with the corresponding optimal
        charge method. It will raise ``RuntimeError`` if the
        FData is not provided by us.
    get_correction()
        Returns `DFT-D3 <https://dftd3.readthedocs.io>`_ correction
        to be used with Fireball. It will return ``None`` if
        either the FData is not provided by us or no optimal parameters are recorded.
    """

    def __init__(self, *, fdata: str, atomsystem: AtomSystem | None = None,
                 fdata_path: str | None = None) -> None:
        self.isloaded = False
        self.wf_loaded = False
        type_check(fdata, str, 'fdata')
        if fdata != 'custom':
            self.name = fdata
            self._from_name()
            self.indb = True
        else:
            type_check(fdata_path, str, 'fdata_path', " when ``fdata='custom'``")
            self.name = fdata
            self.path = fdata_path
            self.indb = False

        if atomsystem is None:
            self.lazy = False
            self.species_present = set()
        else:
            type_check(atomsystem, AtomSystem, 'atomsystem')
            self.lazy = True
            self.species_present = atomsystem.species

        # Fix path
        assert self.path is not None
        if self.path[-1] != os.sep:
            self.path += os.sep
        self.infofile = os.path.join(self.path, 'info.dat')
        if not os.path.isfile(self.infofile):
            raise ValueError("File 'info.dat' not found in the specified "
                            f"fdata path ('{self.path}').")

        # Prepare to (later) create custom info.dat
        with open(self.infofile, 'r') as fp:
            self.infodat = fp.read().splitlines()
        self.blocks = _get_blocks_infodat(self.infodat)
        self._check_species()
        self.pyinfofile = os.path.join(self.path, 'info.py.dat')

    def _from_name(self) -> None:
        self.fdata = get_fdata(self.name)
        self.path = get_data_from_url(self.fdata['url'], self.name, 'FData')
        if self.path[-1] != os.sep:
            self.path += os.sep

    def _check_species(self) -> None:
        self.species = set(_get_specie_block(self.infodat, b) for b in self.blocks)
        if not self.species_present <= self.species:
            raise ValueError(f'FData {self.path} does not contain all required species')

    def _prep_infodat(self) -> None:
        infodat = ['   fireballpy_generated ',
                   f'   {len(self.species) if not self.lazy else len(self.species_present)} - Number of species ']
        for i, j in self.blocks:
            if (not self.lazy) or (_get_specie_block(self.infodat, (i, j)) in self.species_present):
                infodat.extend(self.infodat[i:j+1])
        with open(self.pyinfofile, 'w') as fp:
            fp.write(os.linesep.join(infodat))

    def _find_sizes(self) -> None:
        self.nshells: dict[int, int] = {}
        self.norbitals: dict[int, int] = {}
        with open(self.pyinfofile, 'r') as fp:
            info = fp.read().splitlines()
        n = int(info[1].split('-')[0].strip())
        info = info[2:]  # Remove header
        for _ in range(n):
            z = int(info[3].split('-')[0].strip())
            shellline = info[6]
            shells = list(map(lambda x: int(x.strip()), shellline.split()))
            self.nshells[z] = len(shells)
            self.norbitals[z] = sum([2*l + 1 for l in shells])
            info = info[16:]

    def load_fdata(self) -> None:
        """Set Fortran module variables and load the FData into memory.
        """
        global _loaded_fdata
        assert self.path is not None
        load_tuple = (self.path, self.species_present if self.lazy else self.species)
        if load_tuple != _loaded_fdata:
            self.wf_loaded = False
            self._prep_infodat()
            loadfdata_from_path(self.path)
            _loaded_fdata = deepcopy(load_tuple)
        self._find_sizes()
        self.isloaded = True

    def load_wf(self) -> None:
        if not self.isloaded:
            self.load_fdata()
        if self.wf_loaded:
            return
        assert self.path is not None
        self._prep_infodat()
        self.wfs = {}
        with open(self.pyinfofile, 'r') as fp:
            info = fp.read().splitlines()
        n = int(info[1].split('-')[0].strip())
        info = info[2:]  # Remove header
        for _ in range(n):
            z = int(info[3].split('-')[0].strip())
            self.wfs[z] = {}
            pathline = info[12]
            paths = list(map(lambda x: x.strip(), pathline.split()))
            for path in paths:
                x, y = read_wf(os.path.join(self.path, 'basis', path))
                _, l, _ = read_wf_info(os.path.join(self.path, 'basis', path))
                cs = CubicSpline(x*Bohr, y, extrapolate=False)
                orbname = ANGULAR_MOMENTUM_REV[l]
                if 'e' in path.split('.')[1]:
                    orbname += '*'
                self.wfs[z][orbname] = cs
            info = info[16:]
        self.wf_loaded = True
        os.remove(self.pyinfofile)

    def get_charges_method(self) -> str:
        """Get the optimal charge method for autoconsistency.

        Returns
        -------
        str
            String with the optimal charge method.

        Raises
        ------
        RuntimeError
            If the FData is not provided by us.
        """
        if not self.indb:
            raise RuntimeError('Cannot get charges method for custom FData')
        return self.fdata['charges_method']

    def get_correction(self, charges_method: str) -> dict:
        """Gives the options needed to apply `DFT-D3 <https://dftd3.readthedocs.io>`_ correction.

        Parameters
        ----------
        atomsystem : AtomSystem
            An AtomSystem class with the information of the atomic numbers, positions and the unit cell.

        Returns
        -------
        dict
            Dictionary with the parameters to create the correction.
            Empty if no correction is applicable.
        """
        type_check(charges_method, str, 'charges_method')
        if not self.indb:
            return {}
        try:
            corr = self.fdata['correction']
        except KeyError:
            return {}
        if self.get_charges_method() != charges_method:
            return {}

        return corr


def load_fdata(fdata: str, fdata_path: str | None = None) -> None:
    """Set Fortran module variables and load the whole FData into memory.
    """
    fdatafiles = FDataFiles(fdata=fdata, fdata_path=fdata_path)
    fdatafiles.load_fdata()
