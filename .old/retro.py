from io import StringIO
from os.path import splitext

from ase.data import chemical_symbols  # type: ignore
from ase.io import read  # type: ignore
from ase.io.formats import NameOrFile  # type: ignore


def read_bas(filename: NameOrFile):
    """Reads molecule in .bas format

    Internally this converts the bas file into a IOStream of a xyz file

    Parameters
    ----------
    filename : str or IOStream
        Either the name of the file or a ``IOStream`` with the file.
        Note that the latter is yet not implemented and will raise
        a ``NotImplementedError``.

    Raises
    ------
    NotImplementedError
        If ``filename`` is not a string
    """

    # Not yet implemented load from IO
    if not isinstance(filename, str):
        raise NotImplementedError("Loading from other than a file path \
                is not yet implemented")

    if isinstance(filename, str) and (splitext(filename)[1] != '.bas'):
        raise ValueError(f"{filename} is not a .bas file")

    with open(filename, "r") as fp:
        bas = list(filter(lambda x: x, fp.read().splitlines()))

    if len(bas) != (int(bas[0]) + 1):
        raise ValueError(f"{filename} is not a valid .bas file")

    xyz = [bas[0].strip(), "no comment"]
    for i, line in enumerate(bas[1:]):
        sline = line.split()
        elem = chemical_symbols[int(sline[0])].rjust(3).ljust(8)
        xyz.append(elem + (5*" ").join(f'{float(s):10.5f}' for s in sline[1:]))

    return read(StringIO('\n'.join(xyz)), format='xyz')

# def write_infodat(infodat: InfoDat) -> None:
#    """Writes an InfoDat in the old info.dat format
#
#    Parameters
#    ----------
#    infodat : InfoDat
#        InfoDat to be printed into info.dat file
#    """
#
#    string = []
#    string.append("   pyreball info.dat")
#    string.append(f"  {len(infodat.numbers):2d}\t\t - Number of species")
#    for i, num in enumerate(infodat.numbers):
#        string.append("  " + 70*"=")
#        string.append(f"  {i+1:2d}\t\t - Information for this species")
#        string.append(f"  {infodat.elements[num].symbol:<2}\t\t - Element")
#        string.append(f"  {num:3d}\t\t - Nuclear Z")
#        string.append(
#            f"  {infodat.elements[num].atomic_weight:7.3f}\t - Atomic Mass"
#        )
#        string.append(
#            f"  {infodat.numshells[num]:2d}\t\t - Number of shells; "
#            "L for each shell"
#        )
#        string.append(
#            "  ".join(["", *[f"{x:1d}" for x in infodat.shells[num]]])
#        )
#        string.append(
#            f"  {infodat.numshellsPP[num]:2d}\t\t - Number of shells; "
#            "L for each shell  (Pseudopotential)"
#        )
#        string.append(
#            "  ".join(["", *[f"{x:1d}" for x in infodat.shellsPP[num]]])
#        )
#        string.append(
#            f"  {infodat.cutoffPP[num]:12.4f}\t - Radial cutoffs PP"
#        )
#        string.append(
#            "  ".join(["", *[f"{x:5.2f}" for x in infodat.cutoffs[num]]])
#        )
#        string.append(
#            "  ".join(["", *[f"{x:5.2f}" for x in infodat.qneutral[num]]])
#        )
#        string.append(
#            "  ".join(["", *[f"{x:<25}" for x in infodat.wffiles[num]]])
#        )
#        string.append(
#            "  ".join(["", *[f"{x:<25}" for x in infodat.nafiles[num]]])
#        )
#        string.append(
#            f"  {infodat.energy[num]:12.4f}\t - Atomic energy"
#        )
#        string.append("  " + 70*"=")
#
#        with open("info.dat", "w") as fp:
#            fp.write("\n".join(string))
#
#
# def read_infodat() -> InfoDat:
#    """Loads a legacy info.dat and returns a new InfoDat object
#    """
#
#    with open("info.dat", "r") as fp:
#        dat = fp.readlines()
#
#    nsp = int(dat[1].split("-")[0])
#    shells = {}
#    shellsPP = {}
#    cutoffsPP = {}
#    cutoffs = {}
#    qneutral = {}
#    wffiles = {}
#    nafiles = {}
#    energy = {}
#
#    for i in range(nsp):
#        ii = 16*i + 5
#        n = int(dat[ii].split("-")[0])
#        shells[n] = [
#            int(x) for x in filter(lambda y: y, dat[ii+3][:-1].split(" "))
#        ]
#        shellsPP[n] = [
#            int(x) for x in filter(lambda y: y, dat[ii+5][:-1].split(" "))
#        ]
#        cutoffsPP[n] = float(dat[ii+6].split("-")[0])
#        cutoffs[n] = [
#            float(x) for x in filter(lambda y: y, dat[ii+7][:-1].split(" "))
#        ]
#        qneutral[n] = [
#            float(x) for x in filter(lambda y: y, dat[ii+8][:-1].split(" "))
#        ]
#        wffiles[n] = list(filter(lambda y: y, dat[ii+9][:-1].split(" ")))
#        nafiles[n] = list(filter(lambda y: y, dat[ii+10][:-1].split(" ")))
#        energy[n] = float(dat[ii+11].split("-")[0])
#
#    return InfoDat(shells,
#                   shellsPP,
#                   cutoffsPP,
#                   cutoffs,
#                   qneutral,
#                   wffiles,
#                   nafiles,
#                   energy)
