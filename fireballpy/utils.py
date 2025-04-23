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
                       dst_have: str | None = None) -> tuple[str, time.struct_time]:

    # See if download is needed
    meta = {}
    if dst_have is not None:
        with open(os.path.join(dst_have, 'meta.json'), 'r') as fp:
            meta = json.load(fp)
        if __fb_version__ == meta['VERSION']:
            return dst_have, time.strptime(meta['TIME'], TIMESTR)

    # Make connection
    try:
        r = requests.get(url, stream=True)
    except requests.exceptions.RequestException:
        if dst_have is None:
            raise ConnectionError(f'Could not fetch {filetype}. Please check connection.')
        return dst_have, time.strptime(meta['TIME'], TIMESTR)
    if r.status_code != 200:
        r.raise_for_status()
        raise RuntimeError(f'Fetching {filetype} failed with status code {r.status_code}')

    # Check if update is neeeded
    last_mod = time.strptime(r.headers['Last-Modified'], TIMESTR)
    if dst_have is not None:
        if last_mod == time.strptime(meta['TIME'], TIMESTR):
            # Update version so no more queries are needed
            if __fb_version__ != meta['VERSION']:
                meta['VERSION'] = __fb_version__
                with open(os.path.join(dst_have, 'meta.json'), 'w') as fp:
                    json.dump(meta, fp)
            return dst_have, time.strptime(meta['TIME'], TIMESTR)
        # Ask for update
        global _quiet
        while not _quiet:
            answer = input(f'New version of {filetype} found. Do you want to update? [y|N|q] ')
            if answer.lower() in ['yes', 'y']:
                break
            if answer.lower() in ['', 'no', 'n']:
                return dst_have, time.strptime(meta['TIME'], TIMESTR)
            if answer.lower() in ['quiet', 'q']:
                _quiet = True
                return dst_have, time.strptime(meta['TIME'], TIMESTR)


    # Save to a temp file in dst to avoid permission errors
    # in Windows and cover bad downloads
    # Cover the absurd case in which all posible names are
    # already taken (not going to happen)
    dst_folder = os.path.split(dst)[0]
    if not os.path.isdir(dst_folder):
        os.makedirs(dst_folder)
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

    return dst, last_mod


def extract_tar(tarname: str, tarpath: str, strtime: time.struct_time) -> None:
    tarf = os.path.join(tarpath, tarname)
    tar = tarfile.open(tarf, 'r:xz')
    tar.extractall(tarpath, filter='data')
    if os.path.exists(tarf):
        os.remove(tarf)

    # Write metadata
    meta = {'TIME': time.strftime(TIMESTR, strtime),
            'VERSION': __fb_version__}
    with open(os.path.join(tarf[:-7], 'meta.json'), 'w') as fp:
        json.dump(meta, fp)


def get_data_from_url(url: str, name: str, filetype: str) -> str:
    fb_home = get_fb_home()
    fdata_home = os.path.join(fb_home, 'fdatas')
    if not os.path.isdir(fdata_home):
        os.makedirs(fdata_home)
    if filetype == 'FData':
        dst_have = os.path.join(fb_home, 'fdatas', name)
    elif filetype == 'Pseudopotential':
        dst_have = os.path.join(fb_home, name)
    else:
        raise ValueError("Invalid filetype")

    dst = dst_have + '.tar.xz'
    if not os.path.isdir(dst_have):
        dst_have = ''

    dst, strtime = download_check_tar(url, dst, filetype, dst_have)
    if dst_have and dst == dst_have:
        return dst_have

    tarpath, tarname = os.path.split(dst)
    extract_tar(tarname, tarpath, strtime)
    return dst[:-7]
