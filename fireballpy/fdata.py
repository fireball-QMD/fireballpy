from __future__ import annotations
from copy import deepcopy
from enum import StrEnum
from pathlib import Path
from typing import Any, Iterable

import json
import sys
if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

import numpy as np

from fireballpy.utils import get_data_from_url, read_wf, read_wf_info, ANGULAR_MOMENTUM_REV
from fireballpy.atoms import AtomSystem
from fireballpy.options import ChargesMethod

_loaded_fdata: tuple[str, set[str]] = ("", set(""))
_fdatas: dict[str, Any] | None = None

def _load_db() -> None:
    global _fdatas
    if _fdatas is not None:
        return
    with open(Path(__file__) / "fdata.toml", "rb") as fp:
        _fdatas = tomllib.load(fp)

class FData(StrEnum):
    BIOLOGY = "biology"
    BIOLOGY_EXTENDED = "biology_extended"
    CHON = "CHON"
    HCSI = "HCSi"
    SI = "Si"

def available_fdatas() -> list[str]:
    _load_db()
    global _fdatas
    for i, k in enumerate(FData):
        if "correction" in _fdatas[k]:
            keys[i] += f" ({_fdatas[k]["correction"]["kind"].upper()} available)"
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
    def _get_db_fdata(self, fdata: FData) -> None:
        _load_db()
        global _fdatas
        current = _fdatas[fdata.value]
        self.species: set[str] = set(current["species"])
        self.charges_method = ChargesMethod[current["charges_method"].upper()]
        self.corr: dict[str, Any] = current.get("correction", {})
        self.path = Path(get_data_from_url(current["url"], fdata.value, "FData"))

    def __init__(self, *, fdata: str | Path | FData, atomsystem: AtomSystem | None = None) -> None:
        if isinstance(fdata, FData):
            self._get_db_fdata(fdata)
        else:
            self.path = Path(fdata)
            self.charges_method = None
            self.corr = {}

        with open(self.path / "info.json", "r") as fp:
            self.info: dict[str, Any] = json.load(fp)

        self.species_available = set(self.info["elements"].keys())
        if atomsystem is None:
            self.species_present = self.species_available
        else:
            self.species_present = atomsystem.species
            if self.species_present > self.species_available:
                raise ValueError("provided FData does not have all required species")
            
    def load_fdata(self) -> None:
        """Set Fortran module variables and load the FData into memory.
        """
        global _loaded_fdata
        load_tuple = (self.path, self.species_present)
        if load_tuple != _loaded_fdata:
            load_fdata()

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
