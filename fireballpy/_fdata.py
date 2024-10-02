from __future__ import annotations
import errno
import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
import time
from typing import Any, Iterable, Optional
import uuid

import requests
from tqdm import tqdm

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

ENV_FB_HOME = "FIREBALL_HOME"
ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"

TIMESTR = '%a, %d %b %Y %H:%M:%S %Z'
SHA256SUMS = 'http://fireball.ftmc.uam.es/sha256sums_fireballpy.txt'

_quiet = False


def _load_fdatas() -> dict[str, Any]:
    fdata_db = os.path.join(os.path.dirname(__file__), 'fdata.toml')
    with open(fdata_db, 'rb') as fp:
        return tomllib.load(fp)


def available_fdatas() -> list[str]:
    fdatas = _load_fdatas()
    keys = list(fdatas.keys())
    for i, k in enumerate(keys):
        if 'correction' in fdatas[k]:
            keys[i] += f' ({fdatas[k]["correction"]["type"].upper()} available)'
    return keys


def get_fb_home() -> str:
    return os.path.expanduser(
        os.getenv(ENV_FB_HOME,
                  os.path.join(os.getenv(ENV_XDG_CACHE_HOME,
                                         DEFAULT_CACHE_DIR), "fireball"))
    ) + os.sep


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


# Reference: httlps://github.com/pytorch/pytorch/blob/main/torch/hub.py
def download_file(url: str, dst: str,
                  vhave: Optional[time.struct_time] = None) -> time.struct_time:
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
        while True and not _quiet:
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


def get_fdata_charge_method(name: str) -> str:
    # Check
    fdatas = _load_fdatas()
    if name not in fdatas:
        raise ValueError(f'FData {name} does not exist')
    return fdatas[name]['charges_method']


def get_fdata(name: str, species: set[str]) -> str:
    # Check
    fdatas = _load_fdatas()
    if name not in fdatas:
        raise ValueError(f'FData {name} does not exist')
    if not all([s in fdatas[name]['species'] for s in species]):
        raise ValueError(f'FData {name} does not contain all required species')

    fb_home = get_fb_home()
    os.makedirs(fb_home, exist_ok=True)  # Ensure folder exists
    fdata_folder = os.path.join(fb_home, name) + os.sep

    # Download if not exists
    tarurl = fdatas[name]['url']
    tarpath = os.path.join(fb_home, f'{name}.tar.gz')
    metafile = os.path.join(fdata_folder, 'meta.json')
    if not os.path.isfile(metafile):
        vnew = download_file(tarurl, tarpath)
        meta = {'NAME': name}
        extract_fdata(tarpath, fdata_folder, meta, vnew)
        return fdata_folder

    # Read metadata
    with open(metafile, 'r') as fp:
        meta = json.load(fp)
    vhave = time.strptime(meta['TIME'], TIMESTR)
    vnew = download_file(tarurl, tarpath, vhave)
    if vnew != vhave:
        extract_fdata(tarpath, fdata_folder, meta, vnew)
    return fdata_folder


def _get_blocks_infodat(infodat: list[str]) -> Iterable[tuple[int, int]]:
    eqsplits = [i for i, l in enumerate(infodat) if l == 70*'=']
    return zip(eqsplits[:-1:2], eqsplits[1::2])


def _get_specie_block(infodat: list[str], block: tuple[int, int]) -> str:
    return infodat[block[0]+2].split('-')[0].strip()


def check_fdata_path(fdata_path, species: set[str]) -> str:
    if not isinstance(fdata_path, str):
        raise ValueError("'fdata_path' must be a string.")
    if not os.path.isfile(os.path.join(fdata_path, 'info.dat')):
        raise ValueError("info.dat file not found in the specified "
                         f"fdata path ({fdata_path}).")
    with open(os.path.join(fdata_path, 'info.dat'), 'r') as fp:
        infodat = fp.read().splitlines()
    idatspcs = [_get_specie_block(infodat, b) for b in _get_blocks_infodat(infodat)]
    if not all([s in idatspcs for s in species]):
        raise ValueError(f'FData {fdata_path} does not contain all required species')
    if fdata_path[-1] != os.sep:
        return fdata_path + os.sep
    return fdata_path


def prep_infodat(fdata_path: str, species: Optional[set[str]], lazy: bool) -> set[str]:
    idat = os.path.join(fdata_path, 'info.dat')
    newidat = os.path.join(fdata_path, 'info.ase.dat')
    with open(idat, 'r') as fp:
        infodat = fp.read().splitlines()
    if not lazy:
        idatspcs = [_get_specie_block(infodat, b) for b in _get_blocks_infodat(infodat)]
        with open(newidat, 'w') as fp:
            fp.write(os.linesep.join(infodat))
        return set(idatspcs)

    assert species is not None
    newspecies = []
    newinfodat = []
    for i, j in _get_blocks_infodat(infodat):
        sp = _get_specie_block(infodat, (i, j))
        if sp not in species:
            continue
        newspecies.append(sp)
        newinfodat.extend(infodat[i:j+1])
    newinfodat = ['   fireballpy_generated ',
                  f'   {len(newspecies)} - Number of species '] + newinfodat
    with open(newidat, 'w') as fp:
        fp.write(os.linesep.join(newinfodat))
    return set(newspecies)


def get_correction(name: str) -> tuple[str, dict]:
    fdatas = _load_fdatas()
    if name not in fdatas:
        return '', {}
    if 'correction' not in fdatas[name]:
        return '', {}
    return fdatas[name]['correction']['type'], fdatas[name]['correction']['parameters']
