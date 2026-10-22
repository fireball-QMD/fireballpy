from __future__ import annotations
from copy import deepcopy
import json
import os
import sys
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import numpy as np
from scipy.interpolate import CubicSpline
from ase.units import Bohr

from fireballpy import __version__ as __fb_version__
from fireballpy._errors import type_check
from fireballpy.utils import get_data_from_url, read_wf, read_wf_info, ANGULAR_MOMENTUM_REV, XC_IDS
from fireballpy.atoms import AtomSystem
from fireballpy._fireball import loadfdata_from_path, set_fdata_info

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
def _info_from_infodat(lines: list[str]) -> dict[str, Any]:
    """Parse a legacy 'info.dat' into the same structure 'info.json' has.

    From here on there is a single shape, so nothing downstream has to know
    which of the two the FData came with.
    """
    if 'XC' not in lines[2]:
        raise ValueError("This 'info.dat' predates the libXC header (number of XC "
                         "potentials, their ids and their weights) and cannot be used "
                         "with this version. Regenerate the FData.")
    names = {v: k for k, v in XC_IDS.items()}
    ids = [int(x) for x in lines[3].split()]
    info: dict[str, Any] = {
        'xc': {'names': [names[i] for i in ids],
               'weights': [float(x) for x in lines[4].split()]},
        'species': [],
    }
    blk = lines[5:]
    for _ in range(int(lines[1].split('-')[0].strip())):
        nssh = int(blk[5].split('-')[0].strip())
        nsshPP = int(blk[7].split('-')[0].strip())
        info['species'].append({
            'symbol': blk[2].split('-')[0].strip(),
            'znuc': int(blk[3].split('-')[0].strip()),
            'nshells': nssh,
            'l': [int(x) for x in blk[6].split()][:nssh],
            'qneutral': [float(x) for x in blk[10].split()][:nssh],
            'rcutoff': [float(x) for x in blk[11].split()][:nssh],
            'pseudopotential': {'nshells': nsshPP,
                                'l': [int(x) for x in blk[8].split()][:nsshPP],
                                'rcutoff': float(blk[9].split('-')[0].strip())},
            # Only the old format lists them; 'load_wf' is the only user.
            'wavefunctions': blk[12].split(),
        })
        blk = blk[16:]
    return info


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
        # 'info.json' is what 'fireballpy.basis' writes nowadays, and it is
        # only there when create finished. FDatas published before it fall
        # back to the old 'info.dat'.
        self.jsonfile = os.path.join(self.path, 'info.json')
        self.infofile = os.path.join(self.path, 'info.dat')
        if os.path.isfile(self.jsonfile):
            with open(self.jsonfile, 'r') as fp:
                self.info: dict[str, Any] = json.load(fp)
        elif os.path.isfile(self.infofile):
            with open(self.infofile, 'r') as fp:
                self.info = _info_from_infodat(fp.read().splitlines())
        else:
            raise ValueError("Neither 'info.json' nor 'info.dat' found in the "
                             f"specified fdata path ('{self.path}').")
        self._check_species()

    def _from_name(self) -> None:
        self.fdata = get_fdata(self.name)
        self.path = get_data_from_url(self.fdata['url'], self.name, 'FData')
        if self.path[-1] != os.sep:
            self.path += os.sep

    def _check_species(self) -> None:
        self.species = set(sp['symbol'] for sp in self.info['species'])
        if not self.species_present <= self.species:
            raise ValueError(f'FData {self.path} does not contain all required species')

    def _species(self) -> list[dict[str, Any]]:
        """The species to load: all of them, or only those present when lazy."""
        return [sp for sp in self.info['species']
                if (not self.lazy) or sp['symbol'] in self.species_present]

    def _set_fdata_info(self) -> None:
        """Hand the FData description over to the Fortran side.

        This replaces writing an 'info.py.dat' next to the FData for Fortran
        to parse back: two processes sharing a folder used to race on that
        file, and the layout had to be counted line by line on both sides.
        """
        species = self._species()
        nsh = max(sp['nshells'] for sp in species)
        nspec = len(species)
        znuc = np.array([sp['znuc'] for sp in species], dtype=np.int32)
        nssh = np.array([sp['nshells'] for sp in species], dtype=np.int32)
        nsshpp = np.array([sp['pseudopotential']['nshells'] for sp in species], dtype=np.int32)
        rcpp = np.array([sp['pseudopotential']['rcutoff'] for sp in species], dtype=np.float64)
        lssh = np.zeros((nsh, nspec), dtype=np.int32, order='F')
        lsshpp = np.zeros((nsh, nspec), dtype=np.int32, order='F')
        qneutral = np.zeros((nsh, nspec), dtype=np.float64, order='F')
        rcutoff = np.zeros((nsh, nspec), dtype=np.float64, order='F')
        for i, sp in enumerate(species):
            lssh[:sp['nshells'], i] = sp['l']
            pp = sp['pseudopotential']
            lsshpp[:pp['nshells'], i] = pp['l']
            qneutral[:sp['nshells'], i] = sp['qneutral']
            rcutoff[:sp['nshells'], i] = sp['rcutoff']
        xc = self.info['xc']
        set_fdata_info(znuc, nssh, lssh, nsshpp, lsshpp, rcpp, qneutral, rcutoff,
                       ''.join(f'{n:32s}' for n in xc['names']),
                       np.array(xc['weights'], dtype=np.float64))

    def _find_sizes(self) -> None:
        self.nshells: dict[int, int] = {}
        self.norbitals: dict[int, int] = {}
        for sp in self._species():
            self.nshells[sp['znuc']] = sp['nshells']
            self.norbitals[sp['znuc']] = sum(2*l + 1 for l in sp['l'])

    def load_fdata(self) -> None:
        """Set Fortran module variables and load the FData into memory.
        """
        global _loaded_fdata
        assert self.path is not None
        load_tuple = (self.path, self.species_present if self.lazy else self.species)
        if load_tuple != _loaded_fdata:
            self.wf_loaded = False
            self._set_fdata_info()
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
        self.wfs = {}
        for sp in self._species():
            if 'wavefunctions' not in sp:
                # 'info.json' does not carry the wavefunction filenames:
                # nothing in a calculation reads them, only the plotting in
                # orbitals.py does. Pending rewrite.
                raise NotImplementedError(
                    "load_wf() is not available for FDatas described by 'info.json'; "
                    "it is pending a rewrite.")
            self.wfs[sp['znuc']] = {}
            for path in map(str.strip, sp['wavefunctions']):
                full = os.path.join(self.path, 'basis', path)
                x, y = read_wf(full)
                _, l, _ = read_wf_info(full)
                orbname = ANGULAR_MOMENTUM_REV[l]
                if 'e' in path.split('.')[1]:
                    orbname += '*'
                self.wfs[sp['znuc']][orbname] = CubicSpline(x*Bohr, y, extrapolate=False)
        self.wf_loaded = True

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
