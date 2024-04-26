from __future__ import annotations
from itertools import product
from urllib.request import urlopen
from typing import Iterable

import errno
import os
import shutil
import tempfile
import uuid

from tqdm import tqdm
from fireballpy.infodat import InfoDat

ENV_FB_HOME = "FIREBALL_HOME"
ENV_XDG_CACHE_HOME = "XDG_CACHE_HOME"
DEFAULT_CACHE_DIR = "~/.cache"
ONECFNAMES = ('xc1c_dqi', 'exc1crho', 'nuxc1crho')
TWOCFNAMES = ('overlap', 'vna_ontopl', 'vna_ontopr', 'vna_atom', 'vnl',
              'xc_ontop', 'xc_atom', 'xc_corr', 'dipole_z', 'dipole_y',
              'dipole_x', 'coulomb', 'kinetic', 'nuxc', 'den_ontopl',
              'den_ontopr', 'den_atom', 'dnuxc_ol', 'dnuxc_or',
              'denS_ontopl', 'denS_ontopr', 'denS_atom', 'overlapS')
THREECFNAMES = ('bcna', 'den3', 'deS3')


# Useful wrapper
def _penum(it: Iterable, rep: int = 1) -> product:
    return product(enumerate(it), repeat=rep)


def _get_fb_home() -> str:
    return os.path.expanduser(
        os.getenv(ENV_FB_HOME,
                  os.path.join(os.getenv(ENV_XDG_CACHE_HOME,
                                         DEFAULT_CACHE_DIR), "fireball"))
    ) + os.sep


def _get_downloaded() -> list[str]:
    return [x for x in os.listdir(_get_fb_home())
            if os.path.splitext(x)[1] == ".dat"]


def _get_needed_files(idat: InfoDat) -> list[str]:
    needed = []

    # One centre
    for z in idat.anums:
        needed.extend([f"{x}.{z:02}.dat" for x in ONECFNAMES])

    # Two centres
    for (in1, z1), (in2, z2) in _penum(idat.anums, 2):
        for (i, fname) in enumerate(TWOCFNAMES):
            if i == 13:
                continue

            if i in [1, 14, 17, 19]:
                isub = idat.nshs[z1]
            elif i in [2, 3, 15, 16, 18, 20, 21]:
                isub = idat.nshs[z2]
            elif i in [5, 6, 7]:
                isub = 4
            else:
                isub = 0

            for isorp in range(isub + 1):
                if (i > 13) and (i < 22) and (isorp == 0):
                    continue
                needed.append((f"{fname}_{isorp:02}" if isub else fname) +
                              f".{z1:02}.{z2:02}.dat")

    # Three centres
    for (in1, z1), (in2, z2), (in3, z3) in _penum(idat.anums, 3):
        for isorp, it, (i, fname) in product(range(idat.nshs[z3] + 1),
                                             range(1, 6),
                                             enumerate(THREECFNAMES)):
            if isorp == 0 and i > 0:
                continue
            needed.append(f"{fname}_{it:02}_{isorp:02}"
                          f".{z1:02}.{z2:02}.{z3:02}.dat")

    return needed


# Reference: httlps://github.com/pytorch/pytorch/blob/main/torch/hub.py
def download_file(url: str, dst: str) -> None:
    """Download file from a url and save to a text file in dst

    Parameters
    ----------
    url : str
        URL to the file we want to download
    dst : str
        Path in the hard disk where the file will be saved

    Raises
    ------
    FileExistsError
        If there are no temporary file names available
    """

    u = urlopen(url)
    txt = u.read()
    txt = txt.decode("ascii")

    # Save to a temp file in dst to avoid permission errors
    # in Windows and cover bad downloads
    # Cover the absurd case in which all posible names are
    # already taken (not going to happen)
    dst = os.path.expanduser(dst)
    for _ in range(tempfile.TMP_MAX):
        tmp_dst = f"{dst}.{uuid.uuid4().hex}.tmp"
        try:
            with open(tmp_dst, "x+") as fp:
                fp.write(txt)
        except FileExistsError:
            continue
        break
    else:
        raise FileExistsError(errno.EEXIST, "No temporary file name available")

    # Move to final destination
    try:
        shutil.move(tmp_dst, dst)
    finally:
        if os.path.exists(tmp_dst):
            os.remove(tmp_dst)


def download_needed(idat: InfoDat) -> str:
    """Download FData files which are missing and returns
    their path

    Parameters
    ----------
    idat : InfoDat
        InfoDat object with the information of needed species

    Returns
    -------
    str
        Path to the downloaded FData files
    """
    fb_home = _get_fb_home()
    os.makedirs(fb_home, exist_ok=True)  # Ensure folder exists

    needed = _get_needed_files(idat)
    have = _get_downloaded()

    # Remove old info.dat
    if "info.dat" in have:
        os.remove(os.path.join(fb_home, "info.dat"))

    # Download files we don´t have
    to_download = set(needed) - set(have)
    if len(to_download) > 0:
        for file in tqdm(to_download, unit=" files"):
            download_file(f"https://fireball.ftmc.uam.es/BASE/{file}",
                          os.path.join(fb_home, file))

    # Finally save info.dat
    idat.write_ascii(os.path.join(fb_home, "info.dat"))

    return fb_home


def get_default_infodat() -> InfoDat:
    fb_home = _get_fb_home()
    os.makedirs(fb_home, exist_ok=True)  # Ensure folder exists
    download_file("https://fireball.ftmc.uam.es/BASE/info.dat",
                  os.path.join(fb_home, "default_infodat"))

    return InfoDat(os.path.join(fb_home, "default_infodat"))
