from __future__ import annotations
from copy import deepcopy
import errno
import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
import time
from typing import Any, Iterable
import uuid

import requests
from tqdm import tqdm

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from ._errors import type_check
from ._correction import Correction
from .atoms import AtomSystem
from _fireball import loadfdata_from_path

ENV_FB_HOME = "FIREBALL_HOME"
ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"

TIMESTR = '%a, %d %b %Y %H:%M:%S %Z'
SHA256SUMS = 'https://fireball.ftmc.uam.es/fireballpy/fdata/sha256sums_fireballpy.txt'

_quiet = False
_loaded_fdata = ('', '')


#######################################
# UTIL FUNCTIONS FOR DOWNLOADING FDATAS
def _load_fdatas() -> dict[str, Any]:
    fdata_db = os.path.join(os.path.dirname(__file__), 'fdata.toml')
    with open(fdata_db, 'rb') as fp:
        return tomllib.load(fp)
#######################################

#######################################
# UTIL FUNCTIONS FOR DOWNLOADING FDATAS
# Reference: stackoverflow.com/questions/3431825
def _file_as_blockiter(tarpath: str, blocksize: int = 65536):
    with open(tarpath, 'rb') as fp:
        block = fp.read(blocksize)
        while len(block) > 0:
            yield block
            block = fp.read(blocksize)


# Reference: stackoverflow.com/questions/3431825
def _hash_bytestr_iter(bytesiter, hasher) -> str:
    for block in bytesiter:
        hasher.update(block)
    return hasher.hexdigest()
#######################################

######################################
# UTIL FUNCTIONS FOR CHECKING INFO.DAT
def _get_blocks_infodat(infodat: list[str]) -> Iterable[tuple[int, int]]:
    eqsplits = [i for i, l in enumerate(infodat) if l == 70*'=']
    return list(zip(eqsplits[:-1:2], eqsplits[1::2]))


def _get_specie_block(infodat: list[str], block: tuple[int, int]) -> str:
    return infodat[block[0]+2].split('-')[0].strip()
######################################

def get_fb_home() -> str:
    return os.path.expanduser(
        os.getenv(ENV_FB_HOME,
                  os.path.join(os.getenv(ENV_XDG_CACHE_HOME,
                                         DEFAULT_CACHE_DIR), "fireball"))
    ) + os.sep

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

# Reference: httlps://github.com/pytorch/pytorch/blob/main/torch/hub.py
def download_file(url: str, dst: str,
                  vhave: time.struct_time | None = None) -> time.struct_time:
    # Make connection
    try:
        r = requests.get(url, stream=True)
    except requests.exceptions.RequestException:
        if vhave is None:
            raise ConnectionError('Could not fetch FData. Please check connection.')
        return vhave
    if r.status_code != 200:
        r.raise_for_status()
        raise RuntimeError(f'Fetching FData failed with status code {r.status_code}')

    # See if download is needed
    last_mod = time.strptime(r.headers['Last-Modified'], TIMESTR)
    if last_mod == vhave:
        return last_mod

    # Ask for update
    if vhave is not None and vhave < last_mod:
        global _quiet
        while not _quiet:
            answer = input("New version of FData found. Do you want to update? [y|N|q] ")
            if answer.lower() in ['yes', 'y']:
                break
            if answer.lower() in ['', 'no', 'n']:
                return vhave
            if answer.lower() in ['quiet', 'q']:
                _quiet = True
                return vhave

    # Save to a temp file in dst to avoid permission errors
    # in Windows and cover bad downloads
    # Cover the absurd case in which all posible names are
    # already taken (not going to happen)
    dst = os.path.expanduser(dst)
    file_size = int(r.headers.get('content-length', 0))
    desc = '(Unknown total file size)' if not file_size else ''
    with tqdm.wrapattr(r.raw, 'read', total=file_size, desc=desc) as r_raw:
        for _ in range(tempfile.TMP_MAX):
            tmp_dst = f'{dst}.{uuid.uuid4().hex}.tmp'
            try:
                with open(tmp_dst, "xb") as fp:
                    shutil.copyfileobj(r_raw, fp)
            except FileExistsError:
                continue
            break
        else:
            raise FileExistsError(errno.EEXIST, "No temporary file name available")

    # Get SHA256 sums
    r = requests.get(SHA256SUMS, stream=True)
    for pair in r.content.decode().splitlines():
        sha256sum, shaname = pair.split('  ')
        if shaname == url.split('/')[-1]:
            break
    else:
        raise RuntimeError("sha256sum not found. Something went wrong in connection.")
    mysha256sum = _hash_bytestr_iter(_file_as_blockiter(tmp_dst), hashlib.sha256())
    if sha256sum != mysha256sum:
        raise RuntimeError("sha256sums failed. Try again or check your connection. "
                            "If the error persists please contact us.")

    # Move to final destination
    try:
        shutil.move(tmp_dst, dst)
    finally:
        if os.path.exists(tmp_dst):
            os.remove(tmp_dst)

    return last_mod


def extract_fdata(tarpath: str, fdata_path: str, meta: dict[str, Any],
                  strtime: time.struct_time) -> None:
    tar = tarfile.open(tarpath, 'r:gz')
    tar.extractall(fdata_path, filter='data')
    if os.path.exists(tarpath):
        os.remove(tarpath)

    # Write metadata
    meta['TIME'] = time.strftime(TIMESTR, strtime)
    with open(os.path.join(fdata_path, 'meta.json'), 'w') as fp:
        json.dump(meta, fp)


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
    fdata_path : str | None
        Path to a custom FData. Ignored unless ``fdata = 'custom'``.

    Methods
    -------
    load_fdata(sps, lazy)
        Set Fortran module variables and load the FData
        into memory.
    get_charges_method()
        Return the string with the corresponding optimal
        charge method. It will raise ``RuntimeError`` if the
        FData is not provided by us.
    get_correction(atomsystem)
        Returns `DFT-D3 <https://dftd3.readthedocs.io>`_ correction
        to be used with Fireball. It will return ``None`` if
        either the FData is not provided by us or no optimal parameters are recorded.
    """

    fb_home = get_fb_home()

    def __init__(self, *, fdata: str, fdata_path: str | None = None) -> None:
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

        assert self.path is not None
        if self.path[-1] != os.sep:
            self.path += os.sep
        self.infofile = os.path.join(self.path, 'info.dat')
        if not os.path.isfile(self.infofile):
            raise ValueError("File 'info.dat' not found in the specified "
                            f"fdata path ('{self.path}').")
        with open(self.infofile, 'r') as fp:
            self.infodat = fp.read().splitlines()
        self.blocks = _get_blocks_infodat(self.infodat)
        self.pyinfofile = os.path.join(self.path, 'info.py.dat')

    def _from_name(self) -> None:
        self.fdata = get_fdata(self.name)
        os.makedirs(self.fb_home, exist_ok=True)  # Ensure folder exists
        self.path = os.path.join(self.fb_home, self.name) + os.sep

        # Download if not exists
        tarurl = self.fdata['url']
        tarpath = os.path.join(self.fb_home, f'{self.name}.tar.gz')
        metafile = os.path.join(self.path, 'meta.json')
        if not os.path.isfile(metafile):
            vnew = download_file(tarurl, tarpath)
            meta = {'NAME': self.name}
            extract_fdata(tarpath, self.path, meta, vnew)
            return

        # Read metadata
        with open(metafile, 'r') as fp:
            meta = json.load(fp)
        vhave = time.strptime(meta['TIME'], TIMESTR)
        vnew = download_file(tarurl, tarpath, vhave)
        if vnew != vhave:
            extract_fdata(tarpath, self.path, meta, vnew)

    def _check_species(self, sps: set[str]) -> None:
        self.species = set(_get_specie_block(self.infodat, b) for b in self.blocks)
        if not sps <= self.species:
            raise ValueError(f'FData {self.path} does not contain all required species')

    def _prep_infodat(self, sps: set[str], lazy: bool) -> None:
        self._check_species(sps)
        if lazy:
            self.species = sps
        infodat = []
        for i, j in self.blocks:
            if (not lazy) or (_get_specie_block(self.infodat, (i, j)) in self.species):
                infodat.extend(self.infodat[i:j+1])
        infodat = ['   fireballpy_generated ',
                  f'   {len(self.species)} - Number of species '] + infodat
        with open(self.pyinfofile, 'w') as fp:
            fp.write(os.linesep.join(infodat))

    def load_fdata(self, atomsystem: AtomSystem, lazy: bool) -> None:
        """Set Fortran module variables and load the FData into memory.

        Parameters
        ----------
        atomsystem : AtomSystem
            AtomSystem object with the species to appear in the calculation.
        lazy : bool
            If ``True`` load only FData for the involved species.
            If ``False`` load the full FData.
        """
        global _loaded_fdata
        load_tuple = (self.path, atomsystem.species)
        if load_tuple != _loaded_fdata:
            self._prep_infodat(atomsystem.species, lazy)
            loadfdata_from_path(self.path)
            _loaded_fdata = deepcopy(load_tuple)

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

    def get_correction(self, atomsystem: AtomSystem, charges_method: str) -> tuple[dict, Correction] | tuple[dict, None]:
        """Gives the options needed to apply `DFT-D3 <https://dftd3.readthedocs.io>`_ correction.

        Parameters
        ----------
        atomsystem : AtomSystem
            An AtomSystem class with the information of the atomic numbers, positions and the unit cell.

        Returns
        -------
        DFTD3Correction | None
            Object which can apply the correction in Fireball.
            Will return ``None`` if either we did not provide
            the FData or there are no optimal parameters recorded.
        """
        type_check(atomsystem, AtomSystem, 'atomsystem')
        type_check(charges_method, str, 'charges_method')
        if not self.indb:
            return {}, None
        try:
            corr = self.fdata['correction']
        except KeyError:
            return {}, None
        if self.get_charges_method() != charges_method:
            return {}, None

        return corr, Correction(atomsystem=atomsystem, **corr)
