from __future__ import annotations
import errno
import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import time
import uuid
from typing import Any

import requests
from tqdm import tqdm

from fireballpy import __version__ as __fb_version__


ENV_FB_HOME = "FIREBALL_HOME"
ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"

TIMESTR = '%a, %d %b %Y %H:%M:%S %Z'
SHA256SUMS = 'https://fireball.ftmc.uam.es/fireballpy/sha256sums_fireballpy.txt'

_quiet = False


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


def get_fb_home() -> str:
    return os.path.expanduser(
        os.getenv(ENV_FB_HOME,
                  os.path.join(os.getenv(ENV_XDG_CACHE_HOME,
                                         DEFAULT_CACHE_DIR), "fireball"))
    ) + os.sep


# Reference: https://github.com/pytorch/pytorch/blob/main/torch/hub.py
def download_check_tar(url: str, dst: str, filetype: str,
                       vhave: time.struct_time | None = None) -> time.struct_time:
    # Make connection
    try:
        r = requests.get(url, stream=True)
    except requests.exceptions.RequestException:
        if vhave is None:
            raise ConnectionError(f'Could not fetch {filetype}. Please check connection.')
        return vhave
    if r.status_code != 200:
        r.raise_for_status()
        raise RuntimeError(f'Fetching {filetype} failed with status code {r.status_code}')

    # See if download is needed
    last_mod = time.strptime(r.headers['Last-Modified'], TIMESTR)
    if last_mod == vhave:
        return last_mod

    # Ask for update
    if vhave is not None and vhave < last_mod:
        global _quiet
        while not _quiet:
            answer = input(f'New version of {filetype} found. Do you want to update? [y|N|q] ')
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


def extract_tar(tarpath: str, fdata_path: str, meta: dict[str, Any],
                strtime: time.struct_time) -> None:
    tar = tarfile.open(tarpath, 'r:gz')
    tar.extractall(fdata_path, filter='data')
    if os.path.exists(tarpath):
        os.remove(tarpath)

    # Write metadata
    meta['TIME'] = time.strftime(TIMESTR, strtime)
    meta['VERSION'] = __fb_version__
    with open(os.path.join(fdata_path, 'meta.json'), 'w') as fp:
        json.dump(meta, fp)

