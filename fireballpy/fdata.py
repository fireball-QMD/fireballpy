from __future__ import annotations
from copy import deepcopy
from enum import Enum
from pathlib import Path
from typing import Any

import json
import sys
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import numpy as np

from ase.data import atomic_numbers

from fireballpy import __version__ as __fb_version__
from fireballpy._errors import type_check
from fireballpy.utils import get_data_from_url
from fireballpy.atoms import AtomSystem
from fireballpy.options import ChargesMethod
from fireballpy.basis.fdata import XC_NAMES_LEN
from fireballpy._fireball import loadfdata_from_path, set_fdata_info

_loaded_fdata: tuple[str, list[str]] = ("", [""])
_fdatas: dict[str, Any] = {}

def _load_db() -> None:
    global _fdatas
    if _fdatas:
        return
    with (Path(__file__) / "fdata.toml").open('rb') as fp:
        _fdatas = tomllib.load(fp)

class FData(Enum):
    BIOLOGY = "biology"
    BIOLOGY_EXTENDED = "biology_extended"
    CHON = "CHON"
    HCSI = "HCSi"
    SI = "Si"

def available_fdatas() -> list[str]:
    _load_db()
    keys = []
    for k in FData:
        keys.append(k.value)
        if "correction" in _fdatas[k.value]:
            keys[-1] += f" ({_fdatas[k.value]['correction']['kind'].upper()} available)"
    return keys

class FDataFiles:
    """Class to contain needed information about FData files in disk

    Parameters
    ----------
    fdata : str | pathlib.Path | FData
        Either a provided FData with FData("name"), see ``FData.available_fdatas()``
        for a table of available FData, or the path to a custom FData.
        Please note that the first time a provided FData is used it needs to download
        all the necessary files.
    atomsystem : AtomSystem | None, optional
        AtomSystem object with the information of the species which intervine in the computation.
        If ``None`` (default) then the whole FData will be loaded. This will add a significant overhead
        at the exchange of not having to load new species if they are added later.

    Attributes
    ----------
    name : str
        Name of the provided FData.
    path : pathlib.Path
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
    def _get_db_fdata(self, fdata: FData) -> None:
        _load_db()
        global _fdatas
        current: dict[str, Any] = _fdatas[fdata.value]
        self.charges_method = ChargesMethod[current["charges_method"].upper()]
        self.corr: dict[str, Any] = current.get("correction", {})
        self.path = Path(get_data_from_url(current["url"], fdata.value, "FData"))

    def __init__(self, *, fdata: str | Path | FData, atomsystem: AtomSystem | None = None) -> None:
        if isinstance(fdata, FData):
            self._get_db_fdata(fdata)
        else:
            self.path = Path(fdata)
            self.charges_method = None
            self.corr: dict[str, Any] = {}
        self.path = self.path.absolute()

        with (self.path / "info.json").open('r') as fp:
            self.info: dict[str, Any] = json.load(fp)

        _species: list[tuple[str, int]] = sorted([(sp["symbol"], sp["znuc"]) for sp in self.info["species"]],
                                                 key=lambda x: x[1])
        self.species_available = [x[0] for x in _species]
        if atomsystem is None:
            self.species_present = self.species_available
        else:
            self.species_present = sorted(list(atomsystem.species), key=lambda x: atomic_numbers[x])
            if not set(self.species_present) <= set(self.species_available):
                raise ValueError("provided FData does not have all required species")

    def _set_fdata_info(self) -> None:
        """Hand the FData description over to the Fortran side."""
        species: list[dict[str, Any]] = self.info["species"]
        znuc = np.ascontiguousarray([sp["znuc"] for sp in species], dtype=np.int32)
        nssh = np.ascontiguousarray([len(sp["angular_momentum"]) for sp in species], dtype=np.int32)
        nsshpp = np.ascontiguousarray([len(sp["pseudopotential"]["angular_momentum"]) for sp in species], dtype=np.int32)
        nspec = np.int32(len(species))
        nsh = np.max(nssh)
        lssh = np.zeros((nspec, nsh), dtype=np.int32, order='C')
        rcutoff = np.zeros((nspec, nsh), dtype=np.float64, order='C')
        qneutral = np.zeros((nspec, nsh), dtype=np.float64, order='C')
        qrefxc = np.zeros((nspec, nsh), dtype=np.float64, order='C')
        lsshpp = np.zeros((nspec, nsh), dtype=np.int32, order='C')
        rcpp = np.zeros(nspec, dtype=np.float64, order='C')
        cl_pp = np.zeros((nspec, nsh), dtype=np.float64, order='C')
        for i, sp in enumerate(species):
            rcpp[i] = sp["pseudopotential"]["rcutoff"]
            for j in range(nssh[i]):
                lssh[i, j] = sp["angular_momentum"][j]
                rcutoff[i, j] = sp["rcutoff"][j]
                qneutral[i, j] = sp["qneutral"][j]
                qrefxc[i, j] = sp["qrefxc"][j]
            for j in range(nsshpp[i]):
                lsshpp[i, j] = sp["pseudopotential"]["angular_momentum"][j]
                cl_pp[i, j] = sp["pseudopotential"]["cl"][j]
        exc_names: list[str] = [name.ljust(XC_NAMES_LEN) for name in self.info["xc"]["names"]]
        exc_weights = np.ascontiguousarray(self.info["xc"]["weights"], dtype=np.float64)
        set_fdata_info(znuc, nssh, nsshpp, lssh.T, lsshpp.T, rcutoff.T, qneutral.T, qrefxc.T,
                       rcpp, cl_pp.T, exc_names, exc_weights)
        self.nshells = {int(z): int(n) for z, n in zip(znuc, nssh)}
        self.norbitals = {int(z): int(np.sum([2*lssh[i, j] + 1 for j in range(nssh[i])])) for i, z in enumerate(znuc)}

    def load_fdata(self, lazy: bool = True) -> None:
        """Set Fortran module variables and load the FData into memory.
        """
        if lazy:
            self.species = self.species_present
        else:
            self.species = self.species_available
        global _loaded_fdata
        load_tuple = (str(self.path), self.species)
        if load_tuple != _loaded_fdata:
            self._set_fdata_info()
            loadfdata_from_path(str(self.path / "_")[:-1])
            _loaded_fdata = deepcopy(load_tuple)

    def get_charges_method(self) -> ChargesMethod:
        """Get the optimal charge method for autoconsistency.

        Returns
        -------
        ChargesMethod
            Optimal charge method.

        Raises
        ------
        RuntimeError
            If the FData is not provided by us.
        """
        if self.charges_method is None:
            raise RuntimeError("Cannot get charges method for custom FData")
        return self.charges_method

    def get_correction(self, charges_method: ChargesMethod) -> dict[str, Any]:
        """Gives the options needed to apply `DFT-D3 <https://dftd3.readthedocs.io>`_ correction.

        Parameters
        ----------
        charges_method : ChargesMethod
            Charge projection used.

        Returns
        -------
        dict
            Dictionary with the parameters to create the correction.
            Empty if no correction is applicable.
        """
        if self.charges_method != charges_method:
            return {}
        return self.corr
