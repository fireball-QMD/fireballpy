import json
import os
import pathlib
import shutil
import subprocess
import sys
import re
import time
from enum import IntEnum
from multiprocessing import cpu_count
from typing import Annotated

import numpy as np
from numpy.typing import NDArray

from cyclopts import App, Parameter, Group
from cyclopts.validators import Path
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.text import Text

from fireballpy.utils import TIMESTR, get_fb_home, download_check_tar, extract_tar

from ase.data import atomic_numbers, atomic_names, atomic_masses, chemical_symbols

from fireballpy import __version__ as __fb_version__
from fireballpy.basis import __file__ as __basis_file__
from fireballpy.basis._begin import generate_wavefunctions, generate_vnn


# Commands preparation
console = Console()
app = App(console=console, version=__fb_version__)

## Wavefunctions groups
pot_help = r'The confining potential has the form: V = [r > r_0] * V_0 * exp[-(r_c - r_0)/(r - r_0)/(r_c - r)]'
args_grp = Group.create_ordered("Basic arguments")
parms_grp = Group.create_ordered("Orbital parameters", default_parameter=Parameter(negative='', consume_multiple=True))
exc_grp = Group.create_ordered("Excited states (default: none)")
pot_grp = Group.create_ordered("Confining potential (default: none)", help=pot_help)


def info(msg: str):
    panel = Panel(
        Text(msg, "default"),
        title="Information",
        box=box.ROUNDED,
        expand=True,
        title_align="left",
        style="cyan",
    )
    console.print(panel)

def raise_err(err: Exception):
    panel = Panel(
        Text(str(err), "default"),
        title="Error",
        box=box.ROUNDED,
        expand=True,
        title_align="left",
        style="red",
    )
    console.print(panel)
    exit(1)


PPURL = 'https://fireball.ftmc.uam.es/fireballpy/ppfiles.tar.gz'
ANGULAR_MOMENTUM = {'s': 0, 'p': 1, 'd': 2, 'f': 3}
A0 = {'s': 2.0, 'p': 1.0, 'd': 0.8, 'f': 0.7}
ELECTRONS_PER_ORB = {'s': 2, 'p': 6, 'd': 10, 'f': 14}
ELECTRONIC_CONFIG = '1s,2s,2p,3s,3p,4s,3d,4p,5s,4d,5p,6s,4f,5d,6p,7s,5f,6d,7p,8s'.split(',')


class ECEnum(IntEnum):
    LDA = 3
    BLYP = 9


class ExcitedEnum(IntEnum):
    NONE = 0
    NODE = 1
    DMOL = 2
    MIX = 3


class Element:

    def __init__(self) -> None:
        pass

    def _valence_electrons(self) -> None:
        maxn = 0
        zrem = self._nznuc
        ec = {}
        for orb in ELECTRONIC_CONFIG:
            newn = int(orb[0])
            nele = ELECTRONS_PER_ORB[orb[1]]
            if newn > maxn:
                ec = {}
                maxn = newn
            ec[orb[1]] = (newn, min(zrem, nele))
            if nele < zrem:
                zrem -= nele
                continue
            if newn == maxn:
                ecorbs = list(ec.keys())
                for o in ecorbs:
                    if ec[o][0] < maxn:
                        del ec[o]
            break
        self._valence = [float(ec.get(orb, (0, 0))[1]) for orb in ANGULAR_MOMENTUM]

    @property
    def nznuc(self) -> int:
        return self._nznuc

    @nznuc.setter
    def nznuc(self, value: str):
        self.provided = value
        err = f'Invalid value "{value}" for "ELEMENT" / "--element".'
        if not value.isnumeric():
            if value not in atomic_numbers:
                raise_err(ValueError(f'{err} Unrecognized element.'))
            value = str(atomic_numbers[value])
        self._nznuc = int(value)
        if self._nznuc not in atomic_numbers.values() or self._nznuc < 1:
            raise_err(ValueError(f'{err} Unrecognized element.'))
        self._valence_electrons()

        # No f orbitals
        if self._valence.pop(3) != 0.0:
            raise_err(ValueError(f'{err} Element with unsupported "f" valence orbitals.'))

    @property
    def atomname(self) -> str:
        return atomic_names[self._nznuc]

    @property
    def symbol(self) -> str:
        return chemical_symbols[self._nznuc]

    @property
    def mass(self) -> float:
        return atomic_masses[self._nznuc]

    @property
    def orbitals(self) -> str:
        return self._orbitals

    @orbitals.setter
    def orbitals(self, value: str) -> None:
        if re.match(r'^s?p?d?$', value) is None:
            err = f'Invalid value "{value}" for "ORBITALS" / "--orbitals".'
            raise_err(ValueError(f'{err} Only supported orbitals are "s", "p" and "d". No repetition is supported and they must be provided in order'))
        self._orbitals = value

    @property
    def sav(self) -> NDArray[np.bool_]:
        return self._sav

    @sav.setter
    def sav(self, value: str | None) -> None:
        if value is None:
            self._sav = np.array([True for _ in self._orbitals])
            return
        if re.match(r'^s?p?d?$', value) is None:
            err = f'Invalid value "{value}" for "-s" / "--save".'
            raise_err(ValueError(f'{err} Only supported orbitals are "s", "p" and "d". No repetition is supported and they must be provided in order'))
        self._sav = np.array([x in value for x in self._orbitals])

    @property
    def lam(self) -> list[int]:
        return [ANGULAR_MOMENTUM[o] for o in self._orbitals]

    @property
    def a0(self) -> list[float]:
        return [A0[o] for o in self._orbitals]

    @property
    def xocc(self) -> NDArray[np.float64]:
        return self._xocc

    @xocc.setter
    def xocc(self, value: list[float]) -> None:
        if len(value) != len(self._orbitals):
            err = f'Invalid value "{value}" for "-n" / "--nelectrons".'
            raise_err(ValueError('Must provide one value per orbital.'))
        for v in value:
            err = f'Invalid value "{v}" for "-n" / "--nelectrons".'
            if v < 0.0:
                raise_err(ValueError(f'{err} Must be >= 0.'))
            if v >= 10.0:
                raise_err(ValueError(f'{err} Must be < 10.'))
        self._xocc = np.ascontiguousarray(value, dtype=np.float64)

    @property
    def nzval(self) -> int:
        return int(self._xocc.sum() + 0.1)

    @property
    def rcutoff(self) -> NDArray[np.float64]:
        return self._rcutoff

    @rcutoff.setter
    def rcutoff(self, value: list[float]) -> None:
        if len(value) != len(self._orbitals):
            err = f'Invalid value "{value}" for "-r" / "--radius".'
            raise_err(ValueError('Must provide one value per orbital.'))
        for v in value:
            err = f'Invalid value "{v}" for "-r" / "--radius".'
            if v < 0.0:
                raise_err(ValueError(f'{err} Must be >= 0.'))
            if v >= 100.0:
                raise_err(ValueError(f'{err} Must be < 100.'))
        self._rcutoff = np.ascontiguousarray(value, dtype=np.float64)

    @property
    def rc_max(self) -> float:
        return max(self._rcutoff)

    @property
    def xocc0(self) -> NDArray[np.float64]:
        return self._xocc0

    @xocc0.setter
    def xocc0(self, value: list[float] | None) -> None:
        valence = []
        for l in ANGULAR_MOMENTUM:
            if l in self._orbitals:
                valence.append(self._valence[ANGULAR_MOMENTUM[l]])
        if value is None:
            self._xocc0 = np.ascontiguousarray(valence, dtype=np.float64)
            return
        err = f'Invalid value "{value}" for "-n0" / "--nelectrons-neutral".'
        if len(value) != len(self._orbitals):
            raise_err(ValueError('Must provide one value per orbital.'))

        zval = int(sum(valence) + 0.1)
        if np.abs(sum(value) - zval) > 1e-8:
            raise_err(ValueError(f'{err} They do not add up to valence electrons ({zval}).'))
        for v in value:
            err = f'Invalid value "{v}" for "-n0" / "--nelectrons-neutral".'
            if v < 0.0:
                raise_err(ValueError(f'{err} Must be >= 0.'))
            if v >= 10.0:
                raise_err(ValueError(f'{err} Must be < 10.'))
        self._xocc0 = np.ascontiguousarray(value, dtype=np.float64)

    @property
    def ioption(self) -> int:
        return self._ioption

    @ioption.setter
    def ioption(self, value: int):
        self._ioption = value

    @property
    def nexcite(self) -> int:
        return self._nexcite

    @nexcite.setter
    def nexcite(self, value: int):
        self._nexcite = value

    @property
    def xocc_ion(self) -> NDArray[np.float64]:
        return self._xocc_ion

    @xocc_ion.setter
    def xocc_ion(self, value: list[float] | None):
        if value is None:
            if self._nexcite == 3:
                raise_err(TypeError('If "--excite mix" is specified, then "--nmix" must be provided.'))
            elif self._nexcite == 2:
                if len(self._orbitals) == 1:
                    self._xocc_ion = self.xocc.copy()
                else:
                    self._xocc_ion = self.xocc - 1.0
            else:
                self._xocc_ion = np.zeros(len(self._orbitals), dtype=np.float64, order='C')
            return
        if len(value) != len(self._orbitals):
            err = f'Invalid value "{value}" for "--nmix".'
            raise_err(ValueError(f'{err} Must provide one value per orbital.'))
        for v in value:
            err = f'Invalid value "{v}" for "--nmix".'
            if v < 0.0:
                raise_err(ValueError(f'{err} Must be >= 0.'))
            if v >= 10.0:
                raise_err(ValueError(f'{err} Must be < 10.'))
        self._xocc_ion = np.ascontiguousarray(value, dtype=np.float64)

    @property
    def nzval_ion(self) -> int:
        return int(self._xocc_ion.sum() + 0.1)

    @property
    def cmix(self) -> NDArray[np.float64]:
        return self._cmix

    @cmix.setter
    def cmix(self, value: list[float] | None):
        if value is None:
            if self._nexcite == 3:
                raise_err(TypeError('If "--excite mix" is specified, then "--pmix" must be provided.'))
            self._cmix = np.zeros(len(self._orbitals), dtype=np.float64, order='C')
            return
        if len(value) != len(self._orbitals):
            err = f'Invalid value "{value}" for "--pmix".'
            raise_err(ValueError(f'{err} Must provide one value per orbital.'))
            exit(1)
        for v in value:
            err = f'Invalid value "{v}" for "--pmix".'
            if v < 0.0:
                raise_err(ValueError(f'{err} Must be >= 0.'))
            if v > 1.0:
                raise_err(ValueError(f'{err} Must be <= 1.'))
        self._cmix = np.ascontiguousarray(value, dtype=np.float64)

    def pp_filenames(self, output: pathlib.Path) -> None:
        symchar = str(self._nznuc).rjust(3, '0')

        # Download if needed
        self.ppfile = f'{symchar}.pp'
        ppfolder = os.path.join(get_fb_home(), 'ppfiles')
        pptar = os.path.join(get_fb_home(), 'ppfiles.tar.gz')
        ppfile = os.path.join(ppfolder, 'fbppfiles', self.symbol, str(self.ioption), self.ppfile)
        metafile = os.path.join(ppfolder, 'meta.json')
        if not os.path.isfile(ppfile):
            vnew = download_check_tar(PPURL, pptar, 'PP files')
            meta = {'NAME': 'ppfiles'}
            extract_tar(pptar, ppfolder + os.path.sep, meta, vnew)
        else:
            with open(metafile, 'r') as fp:
                meta = json.load(fp)
            if meta['VERSION'] == __fb_version__:
                return
            vhave = time.strptime(meta['TIME'], TIMESTR)
            vnew = download_check_tar(PPURL, pptar, 'PP files', vhave)
            if vnew != vhave:
                extract_tar(pptar, ppfolder + os.path.sep, meta, vnew)

        # Use normal file if there is no ion
        self.ppionfile = f'{symchar}++.pp'
        ppionfile = os.path.join(ppfolder, 'fbppfiles', self.symbol, str(self.ioption), self.ppionfile)
        if not os.path.isfile(ppionfile):
            ppionfile = ppfile

        # Copy them
        shutil.copy2(ppfile, output / self.ppfile)
        shutil.copy2(ppionfile, output / self.ppionfile)

    def wf_filenames(self) -> None:
        symchar = str(self._nznuc).rjust(3, '0')
        rcchar = [str(self.rc_max).replace('.', '')[:3].ljust(3, '0')]
        for r in self._rcutoff:
            rcchar.append(str(r).replace('.', '')[:3].ljust(3, '0'))
        self.filename_wf, self.filename_ewf = [], []
        for r, o in zip(rcchar[1:], self._orbitals):
            self.filename_wf.append(f'{symchar}_{r}.wf{ANGULAR_MOMENTUM[o] + 1}')
            self.filename_ewf.append(f'{symchar}_{r}.ewf{ANGULAR_MOMENTUM[o] + 1}')

    def na_filenames(self) -> None:
        symchar = str(self._nznuc).rjust(3, '0')
        rcchar = [str(self.rc_max).replace('.', '')[:3].ljust(3, '0')]
        for r in self._rcutoff:
            rcchar.append(str(r).replace('.', '')[:3].ljust(3, '0'))
        self.filename_na, self.filename_ena = [f'{symchar}_{rcchar[0]}.na0'], []
        for r, o in zip(rcchar[1:], self._orbitals):
            self.filename_na.append(f'{symchar}_{r}.na{ANGULAR_MOMENTUM[o] + 1}')
            self.filename_ena.append(f'{symchar}_{r}.ena{ANGULAR_MOMENTUM[o] + 1}')

    def prep_confining(self, vpot: list[float] | None, rpot: list[float] | None) -> None:
        if vpot is None and rpot is None:
            self.v0 = np.zeros(len(self._orbitals), dtype=np.float64, order='C')
            self.r0 = np.zeros(len(self._orbitals), dtype=np.float64, order='C')
            self.ioptim = 0
            return
        if vpot is not None and rpot is not None:
            if len(vpot) != len(self._orbitals):
                err = f'Invalid value "{vpot}" for "--vpot".'
                raise_err(ValueError(f'{err} Must provide one value per orbital.'))
            if len(rpot) != len(self._orbitals):
                err = f'Invalid value "{rpot}" for "--rpot".'
                raise_err(ValueError(f'{err} Must provide one value per orbital.'))
            for v in vpot:
                err = f'Invalid value "{v}" for "--vpot".'
                if v < 0.0:
                    raise_err(ValueError(f'{err} Must be >= 0.'))
            for v in rpot:
                err = f'Invalid value "{v}" for "--rpot".'
                if v < 0.0:
                    raise_err(ValueError(f'{err} Must be >= 0.'))
            self.v0 = np.ascontiguousarray(vpot, dtype=np.float64)
            self.r0 = np.ascontiguousarray(rpot, dtype=np.float64)
            self.ioptim = 1
            return
        if vpot is None:
            raise_err(TypeError('If "--rpot" is specified, then "--vpot" must be provided.'))
        if rpot is None:
            raise_err(TypeError('If "--vpot" is specified, then "--rpot" must be provided.'))


@app.command
def wavefunctions(element: Annotated[str, Parameter(group=args_grp)],
                  orbitals: Annotated[str, Parameter(group=args_grp)], *,
                  nelectrons: Annotated[list[float], Parameter(group=parms_grp, name=['-n', '--nelectrons'])],
                  radius: Annotated[list[float], Parameter(group=parms_grp, name=['-r', '--radius'])],
                  nelectrons_neutral: Annotated[list[float] | None, Parameter(group=parms_grp, name=['-n0', '--nelectrons-neutral'])]=None,
                  exchange_correlation: Annotated[ECEnum, Parameter(group=args_grp, name=['-ec', '--exchange-correlation'])]=ECEnum.BLYP,
                  output: Annotated[pathlib.Path, Parameter(validator=Path(exists=False, file_okay=False, dir_okay=True), group=args_grp,
                                                            name=['-o', '--output'])]=pathlib.Path('cinput'),
                  save: Annotated[str | None, Parameter(group=args_grp, name=['-s', '--save'])]=None,
                  excite: Annotated[ExcitedEnum, Parameter(group=exc_grp, show_default=False)]=ExcitedEnum.NONE,
                  nmix: Annotated[list[float] | None, Parameter(group=exc_grp, negative='')]=None,
                  pmix: Annotated[list[float] | None, Parameter(group=exc_grp, negative='')]=None,
                  vpot: Annotated[list[float] | None, Parameter(group=pot_grp, negative='')]=None,
                  rpot: Annotated[list[float] | None, Parameter(group=pot_grp, negative='')]=None):
    """Application to generate wavefunctions for an ELEMENT with specified ORBITALS.

    Parameters
    ----------
    element: str | int
        Name or atomic number of the desired element.
    orbitals: str
        Orbitals that the element will have ('s', 'p' and 'd').
        No repetition is allowed and they must be provided in order.
        Eg: 's', 'spd', 'sd'.
    nelectrons: list[float]
        Number of electrons for each of the orbitals.
    radius: list[float]
        Cutoff radius for each of the orbitals in atomic units.
    nelectrons_neutral: list[float] | None
        Number of electrons in the ground state for each of the orbitals. [default: valence electrons]
    exchange_correlation: ECEnum
        Exchange-correlation (EC) functional that will be used.
    output: pathlib.Path
        Folder where to save generated files. Will create it if it does not exist.
    save: str
        Which orbitals will be saved. Use only for ion computations. [default: all that enter the computation]
    excite: ExcitedEnum
        If desired, how excited states should be computed: adding a node, using DMOL
        or mixing excited and ground states.
    nmix: list[float]
        [Only with '--excite mix'] Number of electrons for each of the excited orbitals.
    pmix: list[float]
        [Only with '--excite mix'] Percentage of mixing of the ground state for each of the excited orbitals.
    vpot: list[float]
        Parameter V_0 in volts for each of the orbitals.
    rpot: list[float]
        Parameter r_0 in atomic units for each of the orbitals.
    """
    ele = Element()
    ele.nznuc = element
    ele.orbitals = orbitals
    ele.xocc = nelectrons
    ele.rcutoff = radius
    ele.xocc0 = nelectrons_neutral
    ele.ioption = exchange_correlation
    ele.nexcite = excite
    ele.xocc_ion = nmix
    ele.cmix = pmix
    ele.sav = save
    ele.prep_confining(vpot, rpot)

    # Finish early if files exist already
    output = output.absolute()
    os.makedirs(output, exist_ok=True)
    metafile = os.path.join(output, 'meta.json')
    meta = {ele.symbol: {}}
    if os.path.isfile(metafile):
        with open(metafile, 'r') as fp:
            meta = json.load(fp)
        if ele.symbol in meta:
            for o, s in zip(ele.orbitals, ele.sav):
                if not s:
                    continue
                if o in meta[ele.symbol]:
                    raise_err(RuntimeError(f'Files for element "{ele.provided}", orbital "{o}" already exist.'))
            if ele.nexcite != meta[ele.symbol]['nexcite']:
                raise_err(ValueError('Excitation method must be consistent in same element computations.'))
            if ele.ioption != meta[ele.symbol]['ioption']:
                raise_err(ValueError('Exchange-correlation functional must be consistent in same element computations.'))
        else:
            meta[ele.symbol] = {}
    ele.pp_filenames(output)
    ele.wf_filenames()

    meta[ele.symbol]['nexcite'] = ele.nexcite
    meta[ele.symbol]['ioption'] = ele.ioption
    meta[ele.symbol]['ppfile'] = ele.ppfile
    meta[ele.symbol]['ppionfile'] = ele.ppionfile
    for i, o in enumerate(ele.orbitals):
        if not ele.sav[i]:
            continue
        meta[ele.symbol][o] = {'nelectrons': ele.xocc[i],
                               'nelectrons_neutral': ele.xocc0[i],
                               'radius': ele.rcutoff[i],
                               'nmix': ele.xocc_ion[i],
                               'pmix': ele.cmix[i],
                               'vpot': ele.v0[i],
                               'rpot': ele.r0[i],
                               'filename_wf': ele.filename_wf[i],
                               'filename_ewf': ele.filename_ewf[i]}

    generate_wavefunctions(ioption_in=ele.ioption,
                           nexcite_in=ele.nexcite,
                           nznuc_in=ele.nznuc,
                           nzval_in=ele.nzval,
                           nzval_ion_in=ele.nzval_ion,
                           ioptim_in=ele.ioptim,
                           atomname_in=ele.atomname,
                           ppfile_in=ele.ppfile,
                           ppionfile_in=ele.ppionfile,
                           sav_in=ele.sav,
                           outpath_in=str(output) + os.path.sep,
                           lam_in=ele.lam,
                           a0_in=ele.a0,
                           rcutoff_in=ele.rcutoff,
                           xocc_in=ele.xocc,
                           xocc0_in=ele.xocc0,
                           xocc_ion_in=ele.xocc_ion,
                           cmix_in=ele.cmix,
                           r0_in=ele.r0,
                           v0_in=ele.v0,
                           filename_wf_in=ele.filename_wf,
                           filename_ewf_in=ele.filename_ewf)

    with open(metafile, 'w') as fp:
        json.dump(meta, fp)

    if save is not None:
        info(f'Wavefunctions for element "{ele.provided}", orbital{'s' if len(ele.orbitals) > 1 else ''} "{ele.orbitals}" saved.\n'
             f'Please remember to run "fdata wavefunctions finnish {ele.provided} -o {output}" once all desired orbitals have been computed.')
    else:
        compute_vnn(ele, output)
        for i, o in enumerate(ele.orbitals):
            meta[ele.symbol]['filename_na0'] = ele.filename_na[0]
            meta[ele.symbol][o]['filename_na'] = ele.filename_na[i+1]
            meta[ele.symbol][o]['filename_ena'] = ele.filename_ena[i]
        with open(metafile, 'w') as fp:
            json.dump(meta, fp)


@app['wavefunctions'].command
def finnish(element: Annotated[str, Parameter(group=args_grp)], *,
            output: Annotated[pathlib.Path, Parameter(validator=Path(exists=True, file_okay=False, dir_okay=True), group=args_grp,
                                                      name=['-o', '--output'])]=pathlib.Path('cinput')):
    """Generate neutral potentials for the wavefunctions of ELEMENT. Only needed if "-s" / "--save" was specified.

    Parameters
    ----------
    element: str | int
        Name or atomic number of the desired element.
    output: pathlib.Path
        Folder where wavefunctions are located. Neutral potentials will be saved there.
    """
    ele = Element()
    ele.nznuc = element

    output = output.absolute()
    metafile = os.path.join(output, 'meta.json')
    if not os.path.isfile(metafile):
        raise_err(ValueError(f'Invalid value {output} for "-o" / "--output". There is no "meta.json" file.'))
    with open(metafile, 'r') as fp:
        meta = json.load(fp)
    if ele.symbol not in meta:
        raise_err(RuntimeError(f'The "meta.json" does not have information about element "{ele.provided}".'))

    orbitals = [o for o in ANGULAR_MOMENTUM if o in meta[ele.symbol]]
    ele.nexcite = meta[ele.symbol]['nexcite']
    ele.ioption = meta[ele.symbol]['ioption']
    ele.ppfile = meta[ele.symbol]['ppfile']
    ele.ppionfile = meta[ele.symbol]['ppionfile']
    ele.xocc0 = [meta[ele.symbol][o]['nelectrons_neutral'] for o in orbitals]
    ele.rcutoff = [meta[ele.symbol][o]['radius'] for o in orbitals]
    ele.filename_wf = [meta[ele.symbol][o]['filename_wf'] for o in orbitals]
    ele.filename_ewf = [meta[ele.symbol][o]['filename_ewf'] for o in orbitals]

    compute_vnn(ele, output)
    for i, o in enumerate(ele.orbitals):
        meta[ele.symbol]['filename_na0'] = ele.filename_na[0]
        meta[ele.symbol][o]['filename_na'] = ele.filename_na[i+1]
        meta[ele.symbol][o]['filename_ena'] = ele.filename_ena[i]
    with open(metafile, 'w') as fp:
        json.dump(meta, fp)


def compute_vnn(ele: Element, output: pathlib.Path) -> None:
    ele.na_filenames()
    generate_vnn(nexcite_in=ele.nexcite,
                 nznuc_in=ele.nznuc,
                 ppfile_in=ele.ppfile,
                 ppionfile_in=ele.ppionfile,
                 outpath_in=str(output) + os.path.sep,
                 lam_in=ele.lam,
                 rcutoff_in=ele.rcutoff,
                 filename_wf_in=ele.filename_wf,
                 filename_ewf_in=ele.filename_ewf,
                 filename_na_in=ele.filename_na,
                 filename_ena_in=ele.filename_ena)


def write_lines(ele: Element) -> str:
    lines = [ele.symbol,
             ele.nznuc,
             ele.mass,
             ele.ppfile,
             ele.filename_na[0],
             (1 + int(ele.nexcite in [1, 2]))*len(ele.orbitals)]
    for o, n0, r, fw, fn in zip(ele.orbitals, ele.xocc0, ele.rcutoff, ele.filename_wf, ele.filename_na[1:]):
        lines.append(ANGULAR_MOMENTUM[o])
        lines.append(n0)
        lines.append(r)
        lines.append(fw)
        lines.append(fn)
    if ele.nexcite in [1, 2]:
        for o, r, fw, fn in zip(ele.orbitals, ele.rcutoff, ele.filename_ewf, ele.filename_ena):
            lines.append(ANGULAR_MOMENTUM[o])
            lines.append('0')
            lines.append(r)
            lines.append(fw)
            lines.append(fn)
    lines += ['1', '0.5', '0.25 0.25 0.25 0.25 0.25 0.25']

    with open(f'{ele.symbol}.input', 'w') as fp:
        fp.write(os.linesep.join(map(lambda x: str(x), lines)))
        fp.write(os.linesep)
    return f'{ele.symbol}.input'


@app.command
def basis(folder: Annotated[pathlib.Path, Parameter(validator=Path(exists=True, file_okay=False,
                                                                   dir_okay=True))]=pathlib.Path('cinput'), *,
          output: Annotated[pathlib.Path, Parameter(validator=Path(exists=False, file_okay=False, dir_okay=False))]=pathlib.Path('coutput'),
          elements: Annotated[list[str] | None, Parameter(show_default=False, negative=False)]=None,
          verbose: Annotated[bool, Parameter(negative=False)]=False,
          njobs: Annotated[int, Parameter()]=1):
    """Generate the FData for the elements with computed wavefunctions in FOLDER. This process may take some time.

    Parameters
    ----------
    folder: pathlib.Path
        Folder where wavefunctions and neutral potentials are located.
    output: pathlib.Path
        Folder where the FData will be generated.
    elements: list[str] | None
        List of elements (either symbol or atomic number) to consider for the base. [default: all computed in "FOLDER"]
    verbose: bool
        Print which file is being written.
    njobs: int
        Number of processors to use. If set to "-1" all available will be used. To use this option "mpirun" must be available.
    """
    if njobs < 0 and njobs != -1:
        raise_err(ValueError('Invalid value for "--njobs". Must be positive or -1.'))

    folder = folder.absolute()
    metafile = os.path.join(folder, 'meta.json')
    if not os.path.isfile(metafile):
        raise_err(ValueError(f'Invalid value "{folder}" for "FOLDER" / "--folder". There is no "meta.json" file.'))
    with open(metafile, 'r') as fp:
        meta = json.load(fp)

    if elements is None:
        elements = [k for k in meta]
    else:
        for e in elements:
            if e not in meta:
                raise_err(RuntimeError(f'Element "{e}" was not finalized.'))

    els = []
    for e in elements:
        el = Element()
        el.nznuc = e
        for ee in els:
            if el.nznuc == ee.nznuc:
                raise_err(ValueError(f'Invalid value "{elements}" for "ELEMENTS" / "--elements". Contains repeated values.'))
        el.orbitals = ''.join(o for o in ANGULAR_MOMENTUM if o in meta[el.symbol])
        el.xocc0 = [meta[el.symbol][o]['nelectrons_neutral'] for o in el.orbitals]
        el.rcutoff = [meta[el.symbol][o]['radius'] for o in el.orbitals]
        el.filename_wf = [meta[el.symbol][o]['filename_wf'] for o in el.orbitals]
        el.filename_ewf = [meta[el.symbol][o]['filename_ewf'] for o in el.orbitals]
        el.filename_na = [meta[el.symbol]['filename_na0']] + [meta[el.symbol][o]['filename_na'] for o in el.orbitals]
        el.filename_ena = [meta[el.symbol][o]['filename_ena'] for o in el.orbitals]
        el.ppfile = meta[el.symbol]['ppfile']
        el.ppionfile = meta[el.symbol]['ppionfile']
        el.nexcite = meta[el.symbol]['nexcite']
        els.append(el)

    output = output.absolute()
    os.makedirs(output, exist_ok=False)

    nowfolder = os.getcwd()
    os.chdir(output)
    with open('create.input', 'w') as fp:
        fp.write(f'{len(els)}{os.linesep}')
        for el in els:
            shutil.copy2(folder/el.ppfile, el.ppfile)
            shutil.copy2(folder/el.ppionfile, el.ppionfile)
            for f in el.filename_wf + el.filename_ewf + el.filename_na + el.filename_ena:
                if os.path.isfile(folder/f):
                    shutil.copy2(folder/f, f)
            fp.write(write_lines(el))
            fp.write(os.linesep)

    os.makedirs('coutput', exist_ok=False)

    auto = njobs == -1
    if auto:
        njobs = cpu_count()
    try:
        exepath = os.path.join(os.path.split(__basis_file__)[0], 'create.x')
        argv = [str(njobs), f'{1 if verbose else 0}']
        if njobs == 1:
            p = subprocess.Popen([exepath] + argv, stdout=sys.stdout, stderr=sys.stderr)
        else:
            try:
                p = subprocess.Popen(['mpirun', '-np', str(njobs), exepath] + argv,
                                     stdout=sys.stdout, stderr=sys.stderr)
            except FileNotFoundError:
                os.chdir(nowfolder)
                shutil.rmtree(output)
                raise_err(RuntimeError('"mpirun" is not available.'))
        p.communicate()  # type: ignore
        if p.returncode != 0 and auto:  # type: ignore
            p = subprocess.Popen(['mpirun', '-np', str(njobs//2), exepath] + argv,
                                    stdout=sys.stdout, stderr=sys.stderr)
            p.communicate()  # type: ignore
        if p.returncode != 0:  # type: ignore
            os.chdir(nowfolder)
            shutil.rmtree(output)
            raise_err(RuntimeError('Basis creation failed!'))
    except KeyboardInterrupt:
        os.chdir(nowfolder)
        shutil.rmtree(output)
        raise

    os.makedirs('basis')
    shutil.copy2(folder/'meta.json', os.path.join('basis', 'meta.json'))
    shutil.move('create.input', os.path.join('basis', 'create.input'))
    for el in els:
        shutil.move(el.ppfile, os.path.join('basis', el.ppfile))
        shutil.move(el.ppionfile, os.path.join('basis', el.ppionfile))
        shutil.move(f'{el.symbol}.input', os.path.join('basis', f'{el.symbol}.input'))
        for f in el.filename_wf + el.filename_ewf + el.filename_na + el.filename_ena:
            if os.path.isfile(f):
                shutil.move(f, os.path.join('basis', f))

    for f in filter(lambda x: x.endswith('.dat'), os.listdir('coutput')):
        shutil.move(os.path.join('coutput', f), f)
    os.removedirs('coutput')
    os.chdir(nowfolder)


def main():
    app()


if __name__ == '__main__':
    main()
