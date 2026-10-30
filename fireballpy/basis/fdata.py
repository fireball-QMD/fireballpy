from __future__ import annotations
from copy import deepcopy
from enum import IntEnum
from multiprocessing import cpu_count
from typing import Annotated, Any
import json
import os
import pathlib
import shutil
import subprocess
import sys
import re
import warnings
from cyclopts import App, Parameter, Group
from cyclopts.validators import Path
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from ase.data import atomic_numbers, chemical_symbols
import numpy as np
import matplotlib.pyplot as plt

from fireballpy import __version__ as __fb_version__
from fireballpy.utils import get_data_from_url, ANGULAR_MOMENTUM
from fireballpy.wavefunctions import WaveFunction
from fireballpy.basis import __file__ as __basis_file__
from fireballpy.basis._valence import VALENCE_ELECTRONS
from fireballpy.basis._wfs import generate_wavefunctions, generate_vnn

INPUT_VERSION = 1
INFO_VERSION = 1
XC_NAMES_LEN = 32
PPURL = "https://fireball.ftmc.uam.es/fireballpy/ppfiles.tar.xz"
REORBS = re.compile(r"^s?p?d?$")

# Commands preparation
console = Console()
app = App(console=console, version=__fb_version__)

# Wavefunctions groups
args_grp = Group.create_ordered("Basic arguments")
parms_grp = Group.create_ordered("Orbital parameters", default_parameter=Parameter(negative=False, consume_multiple=True))
exc_grp = Group.create_ordered("Exchange-correlation (default: PBE)", default_parameter=Parameter(negative=False, consume_multiple=True))
ex_grp = Group.create_ordered("Excited states (default: none)")
pot_grp = Group.create_ordered("Confining potential (default: none)",
                               help=r"The confining potential has the form: V = [r > r_0] * V_0 * exp[-(r_c - r_0)/(r - r_0)/(r_c - r)]")

def _json_dump(obj: Any, indent: int = 0) -> str:
    """Like json.dumps(indent=2), but a list of numbers or strings stays on one
    line. The file is meant to be read by a person as well as by the code."""
    pad = " "*indent
    if isinstance(obj, dict):
        inner = f",{os.linesep}".join(f"{pad}  {json.dumps(k)}: {_json_dump(v, indent + 2)}"
                            for k, v in obj.items())
        return f"{{{os.linesep}{inner}{os.linesep}{pad}}}"
    if isinstance(obj, list):
        if any(isinstance(x, (dict, list)) for x in obj):
            inner = f",{os.linesep}".join(f"{pad}  {_json_dump(x, indent + 2)}" for x in obj)
            return f"[{os.linesep}{inner}{os.linesep}{pad}]"
        return f"[{', '.join(json.dumps(x) for x in obj)}]"
    return json.dumps(obj)

def info(msg: str) -> None:
    panel = Panel(Text(msg, "default"),
                  title="Info",
                  box=box.ROUNDED,
                  expand=True,
                  title_align="left",
                  style="cyan")
    console.print(panel)

def raise_err(err: Exception) -> None:
    panel = Panel(Text(str(err), "default"),
                  title="Error",
                  box=box.ROUNDED,
                  expand=True,
                  title_align="left",
                  style="red")
    console.print(panel)
    sys.exit(1)

class ExcitedEnum(IntEnum):
    NONE = 0
    NODE = 1
    DMOL = 2
    MIX = 3

class XCEnum(IntEnum):
    LDA = 3
    GGA = 6
    BLYP = 9

class Element:
    def __init__(self, *, element: str, orbitals: str, xocc: list[float], rcutoff: list[float],
                 output: pathlib.Path = pathlib.Path("cinput"), sav: str | None = None, xocc0: list[float] | None = None,
                 nzval_pp: float | None = None,
                 nexcite: ExcitedEnum = ExcitedEnum.NONE, xocc_ion: list[float] | None = None, cmix: list[float] | None = None,
                 exc_names: list[str] | None = None, exc_weights: list[float] | None = None,
                 v0: list[float] | None = None, r0: list[float] | None = None) -> None:
        # nznuc, symbol and valence
        err = f"Invalid value '{element}' for 'ELEMENT' / '--element'."
        if element not in atomic_numbers or element == "X":
            raise_err(ValueError(f"{err} Unrecognized element."))
        self.nznuc = np.int32(atomic_numbers[element])
        self.symbol = element
        self.valence = np.ascontiguousarray(VALENCE_ELECTRONS[int(self.nznuc)], dtype=np.int32)
        if self.valence[3] != 0:
            raise_err(ValueError(f"{err} F-Shell elements are unsupported"))

        # output (existance checks are handled by the CLI)
        self.output = output
        output.mkdir(parents=True, exist_ok=True)
        self.outpath = str(self.output / "_")[:-1]

        # orbitals, save and lam
        if REORBS.match(orbitals) is None:
            err = f"Invalid value '{orbitals}' for 'ORBITALS' / '--orbitals'."
            raise_err(ValueError(f"{err} Only supported orbitals are 's', 'p' and 'd'. "
                                 "No repetition is supported and they must be provided in order"))
        self.orbitals = orbitals
        self.norbs = len(orbitals)
        self.lam = np.ascontiguousarray([ANGULAR_MOMENTUM[o] for o in self.orbitals], dtype=np.int32)
        if sav is None:
            sav = orbitals
        if REORBS.match(sav) is None:
            err = f"Invalid value '{sav}' for '-s' / '--save'."
            raise_err(ValueError(f"{err} Only supported orbitals are 's', 'p' and 'd'. "
                                 "No repetition is supported and they must be provided in order"))
        self.sav = [o in sav for o in self.orbitals]

        # xocc, xocc0, nzval, nzval_pp
        err = f"Invalid value '{xocc}' for '-n' / '--nelectrons'."
        if len(xocc) != self.norbs:
            raise_err(ValueError(f"{err} Must provide one value per orbital."))
        self.xocc = np.ascontiguousarray(xocc, dtype=np.float64)
        if np.any(self.xocc < 0.0):
            raise_err(ValueError(f"{err} All values must be positive."))
        err = f"Invalid value '{xocc0}' for '-n0' / '--nelectrons-neutral'."
        if xocc0 is None:
            xocc0 = [self.valence[int(x)] for x in self.lam]
        if len(xocc0) != self.norbs:
            raise_err(ValueError(f"{err} Must provide one value per orbital."))
        self.xocc0 = np.ascontiguousarray(xocc0, dtype=np.float64)
        if np.any(self.xocc0 < 0.0):
            raise_err(ValueError(f"{err} All values must be positive."))
        if not np.allclose(self.xocc0.sum(), sum(self.valence)):
            raise_err(ValueError(f"{err} They do not add up to the valence electrons."))
        self.nzval = np.int32(self.xocc.sum() + 0.1)
        self.nzval_pp = np.float64(-1.0) if nzval_pp is None else np.float64(nzval_pp)

        # rcutoff
        err = f"Invalid value '{rcutoff}' for '-r' / '--radius'."
        if len(rcutoff) != self.norbs:
            raise_err(ValueError(f"{err} Must provide one value per orbital."))
        self.rcutoff = np.ascontiguousarray(rcutoff, dtype=np.float64)
        if np.any(self.rcutoff < 0.0):
            raise_err(ValueError(f"{err} All values must be positive."))

        # nexcite, xocc_ion, cmix, nzval_ion
        self.exnum = nexcite
        self.excite = nexcite.name
        self.nexcite = np.int32(nexcite.value)
        if xocc_ion is None and (nexcite == ExcitedEnum.DMOL or nexcite == ExcitedEnum.MIX):
            raise_err(TypeError("If '--excite dmol/mix' is specified, then '--nion' must be provided"))
        if xocc_ion is None:
            xocc_ion = [0.0] * self.norbs
        err = f"Invalid value '{xocc_ion}' for '--nion'."
        if len(xocc_ion) != self.norbs:
            raise_err(ValueError(f"{err} Must provide one value per orbital."))
        self.xocc_ion = np.ascontiguousarray(xocc_ion, dtype=np.float64)
        if np.any(self.xocc_ion < 0.0):
            raise_err(ValueError("All values must be positive."))
        if cmix is None and nexcite == ExcitedEnum.MIX:
            raise_err(TypeError("If '--excite mix' is specified, then '--pmix' must be provided"))
        if cmix is None:
            cmix = [1.0] * self.norbs
        err = f"Invalid value '{cmix}' for '--pmix'."
        if len(cmix) != self.norbs:
            raise_err(ValueError(f"{err} Must provide one value per orbital."))
        self.cmix = np.ascontiguousarray(cmix, dtype=np.float64)
        if np.any(self.cmix < 0.0):
            raise_err(ValueError(f"{err} All values must be positive."))
        self.nzval_ion = np.int32(self.xocc_ion.sum() + 0.1)

        # exc_names, exc_weights
        if exc_names is None:
            exc_names = ["gga_x_pbe", "gga_c_pbe"]
            if exc_weights is not None:
                raise_err(TypeError("If '--exc' is not specified, then '--exc-weights' must not be"))
        if exc_weights is None:
            exc_weights = [1.0 for _ in exc_names]
        if len(exc_names) != len(exc_weights):
            raise_err(TypeError("'--exc' and  '--exc-weights' have different lengths"))
        ioption = XCEnum.GGA
        if all(x.startswith("lda") for x in exc_names):
            ioption = XCEnum.LDA
        if not all(x.startswith("gga") for x in exc_names):
            warnings.warn("Some specified XC were LDA and others GGA, we will use GGA pseudopotentials", UserWarning)
        if len(exc_names) == 2 and \
            (exc_names[0] != exc_names[1] and exc_names[0] in ["gga_x_b88", "gga_c_lyp"] and exc_names[1] in ["gga_x_b88", "gga_c_lyp"]):
            ioption = XCEnum.BLYP
            warnings.warn("BLYP XC (and pseudopotentials) is not recommended in Fireball as it is not based in the free gas", UserWarning)
        self.exc_names = [x.ljust(XC_NAMES_LEN) for x in exc_names]
        self.exc_weights = np.ascontiguousarray(exc_weights, dtype=np.float64)
        asort = np.argsort(self.exc_names)
        self.exc_names = [self.exc_names[int(x)] for x in asort]
        self.exc_weights = self.exc_weights[asort]

        # ppfile, ppionfile
        symchar = str(self.nznuc).rjust(3, '0')
        self.ppfile = f"{symchar}.pp"
        self.ppionfile = f"{symchar}++.pp"
        ppfolder = pathlib.Path(get_data_from_url(PPURL, 'ppfiles', 'Pseudopotential'))
        ppfile = ppfolder / self.symbol / str(ioption) / self.ppfile
        ppionfile = ppfolder / self.symbol / str(ioption) / self.ppionfile
        if not ppionfile.is_file():
            ppionfile = ppfile
        shutil.copy2(ppfile, output / self.ppfile)
        shutil.copy2(ppionfile, output / self.ppionfile)

        # v0, r0 and ioptim
        if v0 is None and r0 is not None or v0 is not None and r0 is None:
            raise_err(TypeError("If '--vpot'/'--rpot' is specified, then '--rpot'/'--vpot' must be provided."))
        if v0 is None or r0 is None:
            self.v0 = np.ascontiguousarray([0.0 for _ in self.orbitals], dtype=np.float64)
            self.r0 = np.ascontiguousarray([0.0 for _ in self.orbitals], dtype=np.float64)
            self.ioptim = np.int32(0)
        else:
            err = f"Invalid value '{v0}' for '--vpot'."
            if len(v0) != self.norbs:
                raise_err(ValueError(f"{err} Must provide one value per orbital."))
            self.v0 = np.ascontiguousarray(v0, dtype=np.float64)
            if np.any(self.v0 < 0.0):
                raise_err(ValueError(f"{err} All values must be positive."))
            err = f"Invalid value '{r0}' for '--rpot'."
            if len(r0) != self.norbs:
                raise_err(ValueError(f"{err} Must provide one value per orbital."))
            self.r0 = np.ascontiguousarray(r0, dtype=np.float64)
            if np.any(self.r0 < 0.0):
                raise_err(ValueError(f"{err} All values must be positive."))
            self.ioptim = np.int32(1)

    def serialize_and_check(self, metadata: dict[str, Any]) -> dict[str, Any]:
        meta: dict[str, Any] = deepcopy(metadata)
        meta["commands"] = meta.get("commands", [])
        meta["commands"].append(" ".join(sys.argv))
        meta["excite"] = meta.get("excite", self.excite)
        if meta["excite"] != self.excite:
            raise_err(RuntimeError("Cannot use different methods for excitation for the same element"))
        meta["exc"] = meta.get("exc", self.exc_names)
        meta["exc_weights"] = meta.get("exc_weights", self.exc_weights)
        if not np.all(meta["exc"] == self.exc_names) or not np.all(meta["exc_weights"] == self.exc_weights):
            raise_err(RuntimeError("Cannot use different XC for the same element"))
        meta["exc"] = [x.strip() for x in self.exc_names]
        meta["exc_weights"] = [x for x in self.exc_weights]
        meta["ppfile"], meta["ppionfile"] = self.ppfile, self.ppionfile

        meta["orbitals"] = meta.get("orbitals", "")
        for i, o in enumerate(self.orbitals):
            if not self.sav[i]:
                continue
            meta["orbitals"] += o
            meta[o] = {"nelectrons": self.xocc[i], "nelectrons_neutral": self.xocc0[i], "nzval_pp": self.nzval_pp, "radius": self.rcutoff[i]}
            if self.nexcite == ExcitedEnum.DMOL or self.nexcite == ExcitedEnum.MIX:
                meta[o]["nion"] = self.xocc_ion[i]
                if self.nexcite == ExcitedEnum.MIX:
                    meta[o]["pmix"] = self.cmix[i]
        return meta

    @classmethod
    def load(cls, element: str | list[str] | None, opt: str, folder: pathlib.Path = pathlib.Path("cinput")) -> list[Element]:
        metafile = folder / "meta.json"
        if not metafile.is_file():
            raise_err(ValueError(f"Invalid value '{folder}' for '{opt}'. There is no 'meta.json' file."))
        with metafile.open('r') as fp:
            meta: dict[str, dict[str, Any]] = json.load(fp)
        if element is None:
            elements = [k for k in meta]
        elif isinstance(element, str):
            elements = element
        else:
            elements = element
        ret: list[Element] = []
        for el in elements:
            if el not in meta:
                raise_err(RuntimeError(f"The 'meta.json' does not have information about el '{el}'."))
            elmeta = meta[el]
            orbitals = ""
            if "s" in elmeta["orbitals"]:
                orbitals += "s"
            if "p" in elmeta["orbitals"]:
                orbitals += "p"
            if "d" in elmeta["orbitals"]:
                orbitals += "d"
            assert len(orbitals) > 0
            xocc = [float(elmeta[o]["nelectrons"]) for o in orbitals]
            rcutoff = [float(elmeta[o]["radius"]) for o in orbitals]
            xocc0 = [float(elmeta[o]["nelectrons_neutral"]) for o in orbitals]
            xocc_ion = None
            cmix = None
            if elmeta["excite"] == ExcitedEnum.DMOL.name or elmeta["excite"] == ExcitedEnum.MIX.name:
                xocc_ion = [float(elmeta[o]["nion"]) for o in orbitals]
                if elmeta["excite"] == ExcitedEnum.MIX.name:
                    cmix = [float(elmeta[o]["pmix"]) for o in orbitals]
            ret.append(cls(element=el, orbitals=orbitals, output=folder, xocc=xocc, rcutoff=rcutoff, xocc0=xocc0,
                       nzval_pp=elmeta[orbitals[0]]["nzval_pp"], nexcite=ExcitedEnum[elmeta["excite"]],
                       xocc_ion=xocc_ion, cmix=cmix, exc_names=elmeta["exc"], exc_weights=elmeta["exc_weights"]))
        return ret

@app.command
def wavefunctions(element: Annotated[str, Parameter(group=args_grp)],
                  orbitals: Annotated[str, Parameter(group=args_grp)], *,
                  nelectrons: Annotated[list[float], Parameter(group=parms_grp, name=["-n", "--nelectrons"])],
                  radius: Annotated[list[float], Parameter(group=parms_grp, name=["-r", "--radius"])],
                  nelectrons_neutral: Annotated[list[float] | None, Parameter(group=parms_grp, name=["-n0", "--nelectrons-neutral"],
                                                                              show_default="valence electrons")]=None,
                  valence_pp: Annotated[float | None, Parameter(group=parms_grp, name=["-zpp", "--valence-pp"],
                                                                show_default="nuclear number")]=None,
                  output: Annotated[pathlib.Path, Parameter(validator=Path(exists=False, file_okay=False, dir_okay=True), group=args_grp,
                                                            name=["-o", "--output"])]=pathlib.Path("cinput"),
                  save: Annotated[str | None, Parameter(group=args_grp, name=["-s", "--save"])]=None,
                  exc: Annotated[list[str] | None, Parameter(group=exc_grp, show_default="gga_x_pbe gga_c_pbe")]=None,
                  exc_weights: Annotated[list[float] | None, Parameter(group=exc_grp, show_default="1.0 each")]=None,
                  excite: Annotated[ExcitedEnum, Parameter(group=ex_grp, show_default=False)]=ExcitedEnum.NONE,
                  nion: Annotated[list[float] | None, Parameter(group=ex_grp, negative=False, consume_multiple=True)]=None,
                  pmix: Annotated[list[float] | None, Parameter(group=ex_grp, negative=False, consume_multiple=True)]=None,
                  vpot: Annotated[list[float] | None, Parameter(group=pot_grp, negative=False, consume_multiple=True)]=None,
                  rpot: Annotated[list[float] | None, Parameter(group=pot_grp, negative=False, consume_multiple=True)]=None) -> None:
    """Application to generate wavefunctions for an ELEMENT with specified ORBITALS.

    Parameters
    ----------
    element: str
        Name of the desired element.
    orbitals: str
        Orbitals that the element will have ('s', 'p' and 'd').
        No repetition is allowed and they must be provided in order.
        Eg: 's', 'spd', 'sd'.
    nelectrons: list[float]
        Number of electrons for each of the orbitals.
    radius: list[float]
        Cutoff radius for each of the orbitals in atomic units.
    nelectrons_neutral: list[float] | None
        Number of electrons in the ground state for each of the orbitals.
    valence_pp: float | None
        Number of electrons considered for the pseudopotential. Increasing will lead to more constrained orbitals.
    output: pathlib.Path
        Folder where to save generated files. Will create it if it does not exist.
    save: str
        Which orbitals will be saved. Use only for ion computations.
    exc: list[str] | None
        List of exchange-correlation (EXC) functional that will be used following LibXC naming.
    exc_weights: list[float] | None
        Weights associated to each exchange-correlation functional.
    excite: ExcitedEnum
        If desired, how excited states should be computed: adding a node, using DMOL
        or mixing excited and ground states.
    nion: list[float]
        [Required for '--excite dmol/mix'] Number of electrons for each of the excited orbitals.
    pmix: list[float]
        [Required for '--excite mix'] Percentage of mixing of the ground state for each of the excited orbitals.
    vpot: list[float]
        Parameter V_0 in volts for each of the orbitals.
    rpot: list[float]
        Parameter r_0 in atomic units for each of the orbitals.
    """
    el = Element(element=element, orbitals=orbitals, xocc=nelectrons, rcutoff=radius, output=output, sav=save,
                 xocc0=nelectrons_neutral, nzval_pp=valence_pp, nexcite=excite, xocc_ion=nion, cmix=pmix,
                 exc_names=exc, exc_weights=exc_weights, v0=vpot, r0=rpot)

    # Finish early if files exist already
    for i, o in enumerate(el.orbitals):
        if not el.sav[i]:
            continue
        if (output / f"{el.nznuc:03}_{o}.wf").is_file():
            raise_err(RuntimeError(f"Files for element '{el.symbol}', orbital '{o}' already exist"))

    # Prepare serialization and handle inconsistencies
    metafile = output / "meta.json"
    meta: dict[str, dict[str, Any]] = {}
    if metafile.is_file():
        with metafile.open('r') as fp:
            meta = json.load(fp)
    if el.symbol in meta and save is None:
        raise_err(RuntimeError(f"Files for element '{el.symbol}', already exist, please specify '-s'/'--save' explicitely"))
    elmeta: dict[str, Any] = meta.get(el.symbol, {})
    meta[el.symbol] = elmeta | el.serialize_and_check(elmeta)
    generate_wavefunctions(nznuc_in=el.nznuc, nexcite_in=el.nexcite, nzval_in=el.nzval,
                           nzval_ion_in=el.nzval_ion, nzval_pp_in=el.nzval_pp, ioptim_in=el.ioptim,
                           ppfile_in=el.ppfile, ppionfile_in=el.ppionfile, outpath_in=el.outpath,
                           sav_in=el.sav, lam_in=el.lam, rcutoff_in=el.rcutoff,
                           xocc_in=el.xocc, xocc0_in=el.xocc0, xocc_ion_in=el.xocc_ion, cmix_in=el.cmix,
                           exc_names_in=el.exc_names, exc_weights_in=el.exc_weights, r0_in=el.r0, v0_in=el.v0)

    with metafile.open('w') as fp:
        json.dump(meta, fp)
    if save is not None:
        info(f"Wavefunctions for element '{el.symbol}', orbitals '{''.join(o for i, o in enumerate(el.orbitals) if save[i])}' saved.\n"
             f"Please remember to run 'fdata wavefunctions finnish {el.symbol} -o {output}' once all desired orbitals have been computed.")
    else:
        generate_vnn(nexcite_in=el.nexcite, nznuc_in=el.nznuc, ppfile_in=el.ppfile,
                     ppionfile_in=el.ppionfile, outpath_in=el.outpath,
                     lam_in=el.lam, rcutoff_in=el.rcutoff)

@app['wavefunctions'].command
def finnish(element: Annotated[str, Parameter(group=args_grp)], *,
            output: Annotated[pathlib.Path, Parameter(validator=Path(exists=True, file_okay=False, dir_okay=True), group=args_grp,
                                                      name=["-o", "--output"])]=pathlib.Path("cinput")) -> None:
    """Generate neutral potentials for the wavefunctions of ELEMENT. Only needed if "-s" / "--save" was specified.

    Parameters
    ----------
    element: str | int
        Name or atomic number of the desired element.
    output: pathlib.Path
        Folder where wavefunctions are located. Neutral potentials will be saved there.
    """
    el = Element.load(element, opt="'-o' / '--output'", folder=output)[0]
    generate_vnn(nexcite_in=el.nexcite, nznuc_in=el.nznuc, ppfile_in=el.ppfile,
                 ppionfile_in=el.ppionfile, outpath_in=el.outpath,
                 lam_in=el.lam, rcutoff_in=el.rcutoff)

def write_create_input(els: list[Element], path: pathlib.Path, infolder: pathlib.Path,
                       refchs: list[list[float]], exc_names: list[str], exc_weights: list[float]) -> None:
    """Write the single file that ``create`` parses."""
    lines = ["# create.ini - written by fireballpy.basis, to be parsed by create",
             "[meta]",
             f"version   = {INPUT_VERSION}",
             f"generator = fireballpy {__fb_version__}",
             f"folder    = {str(infolder / '_')[:-1]}",
             "",
             "[xc]",
             f"names   = {' '.join(exc_names)}",
             f"weights = {' '.join(f'{w:g}' for w in exc_weights)}",
             "",
             "[species]",
             f"count = {len(els)}"]
    for i, el in enumerate(els):
        shells = [(0, o, float(n0), qref) for o, n0, qref in zip(el.orbitals, el.xocc0, refchs[i][:el.norbs])]
        if el.exnum == ExcitedEnum.NODE or el.exnum == ExcitedEnum.DMOL:
            shells += [(1, o, 0.0, qref) for o, qref in zip(el.orbitals, refchs[i][el.norbs:])]
        lines += ["",
                  f"[species.{i + 1}]",
                  f"symbol  = {el.symbol}",
                  f"znuc    = {el.nznuc}",
                  f"ppfile  = {el.ppfile}",
                  f"nshells = {len(shells)}",
                  "# shell = excited angular_momentum qneutral qrefxc"]
        lines += [f"shell = {ex} {ANGULAR_MOMENTUM[l]} {n0:.6f} {qref:.6f}" for ex, l, n0, qref in shells]
    lines += ""
    path.write_text(os.linesep.join(lines))

def read_pp_data(path: pathlib.Path) -> dict[str, Any]:
    """Shells, angular momenta and cutoff of a fhi98PP pseudopotential."""
    with path.open('r') as fp:
        lines = fp.read().splitlines()
    nshells = int(lines[15].split()[0])
    angular_momentum = [int(x) for x in lines[16].split()[:nshells]]
    rcutoff = float(lines[19].split()[0])
    np = int(lines[20].strip())
    rest = lines[(21 + np):]
    cl: list[float] = []
    for _ in range(nshells):
        np = int(rest[0].split()[2].strip())
        rest = rest[(np + 1):]
    for _ in range(nshells):
        np = int(rest[0].split()[2].strip())
        cl.append(float(rest[0].split()[4].strip()))
        rest = rest[(np + 1):]
    return {"angular_momentum": angular_momentum, "rcutoff": rcutoff, "cl": cl}

def write_info_json(els: list[Element], inpath: pathlib.Path, outpath: pathlib.Path,
                    refchs: list[list[float]], exc_names: list[str], exc_weights: list[float]) -> None:
    """Write what ``fireball`` needs in order to read this FData."""
    info: dict[str, Any] = {
        "version": INFO_VERSION,
        "generator": f"fireballpy {__fb_version__}",
        "xc": {"names": [x.strip() for x in exc_names],
               "weights": [x for x in exc_weights]},
        "species": [],
    }
    for i, el in enumerate(els):
        shells = [(o, float(q), qref, float(r))
                  for o, q, qref, r in zip(el.orbitals, el.xocc0, refchs[i][:el.norbs], el.rcutoff)]
        if el.exnum == ExcitedEnum.NODE or el.exnum == ExcitedEnum.DMOL:
            shells += [(o, 0.0, qref, float(r)) for o, qref, r in zip(el.orbitals, refchs[i][el.norbs:], el.rcutoff)]
        pp = read_pp_data(inpath / el.ppfile)
        info["species"].append({
            "symbol": el.symbol,
            "znuc": int(el.nznuc),
            "angular_momentum": [ANGULAR_MOMENTUM[o] for o, _, _, _ in shells],
            "qneutral": [q for _, q, _, _ in shells],
            "qrefxc": [qref for _, _, qref, _ in shells],
            "rcutoff": [r for _, _, _, r in shells],
            "pseudopotential": {"angular_momentum": pp["angular_momentum"],
                                "rcutoff": pp["rcutoff"],
                                "cl": pp["cl"]},
        })
    with (outpath / "info.json").open('w') as fp:
        fp.write(_json_dump(info))
        fp.write(os.linesep)

@app.command
def basis(folder: Annotated[pathlib.Path, Parameter(validator=Path(exists=True, file_okay=False, dir_okay=True))]=pathlib.Path("cinput"), *,
          output: Annotated[pathlib.Path, Parameter(validator=Path(exists=False, file_okay=False, dir_okay=True))]=pathlib.Path("coutput"),
          elements: Annotated[list[str] | None, Parameter(show_default="all inside 'FOLDER'", negative=False, consume_multiple=True)]=None,
          exc: Annotated[list[str] | None, Parameter(consume_multiple=True, negative=False, show_default="gga_x_pbe gga_c_pbe")]=None,
          exc_weights: Annotated[list[float] | None, Parameter(consume_multiple=True, negative=False, show_default="1.0 each")]=None,
          exc_charges: Annotated[list[str] | None, Parameter(consume_multiple=True, negative=False, show_default="neutral charges")]=None,
          verbose: Annotated[bool, Parameter(negative=False)]=False,
          njobs: Annotated[int, Parameter()]=1) -> None:
    """Generate the FData for the elements with computed wavefunctions in FOLDER. This process may take some time.

    Parameters
    ----------
    folder: pathlib.Path
        Folder where wavefunctions and neutral potentials are located.
    output: pathlib.Path
        Folder where the FData will be generated.
    elements: list[str] | None
        List of elements (either symbol or atomic number) to consider for the base.
    exc: list[str] | None
        List of exchange-correlation (EXC) functional that will be used following LibXC naming.
    exc_weights: list[float] | None
        Weights associated to each exchange-correlation functional.
    exc_charges: list[str] | None
        Reference charges for XC in the format of ELEMENT_1 SHELL_CHARGE_1 SHELL_CHARGE_2 ... ELEMENT_2 ...
    verbose: bool
        Print which file is being written.
    njobs: int
        Number of processors to use. If set to "-1" all available will be used. To use this option "mpirun" must be available.
    """
    if njobs < 0 and njobs != -1:
        raise_err(ValueError("Invalid value for '--njobs'. Must be positive or -1."))
    folder = folder.absolute()
    output = output.absolute()
    els = sorted(Element.load(elements, opt="'FOLDER' / '--folder'", folder=folder), key=lambda el: el.nznuc)
    if exc is None:
        exc = ["gga_x_pbe", "gga_c_pbe"]
        if exc_weights is not None:
            raise_err(TypeError("If '--exc' is not specified, then '--exc-weights' must not be"))
    if exc_weights is None:
        exc_weights = [1.0 for _ in exc]
    if len(exc) != len(exc_weights):
        raise_err(TypeError("'--exc' and  '--exc-weights' have different lengths"))
    if exc_charges is None:
        refchs = [[float(el.xocc0[i]) for i in range(el.norbs)] for el in els]
        for i, el in enumerate(els):
            if el.exnum == ExcitedEnum.NODE or el.exnum == ExcitedEnum.DMOL:
                refchs[i] += [0.0 for _ in range(el.norbs)]
    else:
        split_charges = list(x.strip() for x in exc_charges)
        refchs: list[list[float]] = []
        for el in els:
            if el.symbol not in split_charges:
                refchs.append([float(el.xocc0[i]) for i in range(el.norbs)])
                if el.exnum == ExcitedEnum.NODE or el.exnum == ExcitedEnum.DMOL:
                    refchs[-1] += [0.0 for _ in range(el.norbs)]
                continue
            i = split_charges.index(el.symbol)
            try:
                refchs.append([float(split_charges[i + j + 1]) for j in range(el.norbs)])
                if el.exnum == ExcitedEnum.NODE or el.exnum == ExcitedEnum.DMOL:
                    refchs[-1] += [float(split_charges[i + el.norbs + j + 1]) for j in range(el.norbs)]
            except ValueError:
                raise_err(ValueError("Invalid value for '--exc-charges'. Could not be parsed."))
    output.mkdir(parents=True, exist_ok=True)
    create_input = output / "create.ini"
    write_create_input(els, create_input, folder, refchs, exc, exc_weights)

    njobs = cpu_count() if njobs == -1 else njobs
    try:
        exepath = pathlib.Path(__basis_file__).parent / "create.x"
        if njobs == 1:
            command = [str(exepath)]
            if verbose:
                command += ["-v"]
            p = subprocess.Popen(command, cwd=output, stdout=sys.stdout, stderr=sys.stderr)
        else:
            command = ["mpirun", "-np", f"{njobs}", str(exepath)]
            if verbose:
                command += ["-v"]
            p = subprocess.Popen(command, cwd=output, stdout=sys.stdout, stderr=sys.stderr)
        p.communicate()
        if p.returncode != 0:
            shutil.rmtree(output)
            raise_err(RuntimeError("Basis creation failed!"))
    except (KeyboardInterrupt, FileNotFoundError):
        shutil.rmtree(output)
        raise

    # Only now that create has finished without error.
    write_info_json(els, folder, output, refchs, exc, exc_weights)
    basis_folder = output / "basis"
    basis_folder.mkdir(exist_ok=True)
    for el in els:
        shutil.copy2(folder / el.ppfile, basis_folder)
        shutil.copy2(folder / el.ppionfile, basis_folder)
        for o in el.orbitals:
            shutil.copy2(folder / f"{int(el.nznuc):03}_{o}.wf", basis_folder)
            shutil.copy2(folder / f"{int(el.nznuc):03}_{o}.na", basis_folder)
            shutil.copy2(folder / f"{int(el.nznuc):03}.na0", basis_folder)
            if el.exnum == ExcitedEnum.NODE or el.exnum == ExcitedEnum.DMOL:
                shutil.copy2(folder / f"{int(el.nznuc):03}_{o}.ewf", basis_folder)
                shutil.copy2(folder / f"{int(el.nznuc):03}_{o}.ena", basis_folder)
    shutil.copy2(folder / "meta.json", basis_folder)
    shutil.copy2(create_input, basis_folder)

@app.command
def plot(wavefunctions: Annotated[list[pathlib.Path],
                                  Parameter(negative=False, validator=Path(exists=True, file_okay=True, dir_okay=False),
                                            consume_multiple=True)]):
    """Plot wavefunctions.

    It will generate a matplotlib.pyplot plot of the specified wavefunctions.
    Options for the style of the plot and file outputs are still under development.

    Parameters
    ----------
    wavefunctions: list[pathlib.Path]
        List of wavefunction files (.wf or .ewf) to be plotted. It will raise an exception
        if unknown file extensions are provided.
    """
    excited = []
    for path in wavefunctions:
        if path.suffix == ".ewf":
            excited.append(True)
        elif path.suffix == ".wf":
            excited.append(False)
        else:
            raise_err(ValueError(f"Invalid value '{path}' for 'WAVEFUNCTIONS'. Unrecognized file extension."))

    colors = list(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    my_colors = colors * (len(wavefunctions) // len(colors)) + colors[:len(wavefunctions) % len(colors)]
    _, ax = plt.subplots(figsize=(6.9, 4.3), layout="constrained")
    for i, path in enumerate(wavefunctions):
        wf = WaveFunction.load(path)
        symbol = chemical_symbols[wf.z]
        orbital = wf.lsym

        sub = f"{orbital}{'*' if excited[i] else ''}"
        label = r"$\phi_{" + sub + r"}(" + symbol + r")\ \ R_c=" + f"{wf.rcut:3f}" + r"$"
        ax.plot(wf.r, wf.psi, lw=2, ls='-', color=my_colors[i], label=label)

    ax.axhline(lw=1, color='k')
    ax.set_xlabel("Distance ($a_0$)", fontsize=16)
    ax.set_ylabel(r"$| \phi |$", fontsize=16)
    ax.set_xlim(left=0)
    ax.legend(fontsize=14)
    plt.show()

def main():
    app()

if __name__ == "__main__":
    app()
