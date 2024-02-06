from io import StringIO
from os.path import splitext

from ase.io import read
from ase.io.formats import NameOrFile

from mendeleev import element

def read_bas(filename: NameOrFile):

    # Not yet implemented load from IO
    if not isinstance(filename, str):
        raise NotImplementedError("Loading from other than a file path is not yet implemented")

    if isinstance(filename, str) and (splitext(filename)[1] != '.bas'):
        raise ValueError(f"{filename} is not a .bas file")

    with open(filename, "r") as f:
        bas = f.readlines()

    while bas[-1] == '\n':
        bas = bas[:-1]

    if len(bas) != (int(bas[0]) + 1):
        raise ValueError(f"{filename} is not a valid .bas file")
    xyz = [bas[0].rstrip(), "no comment"]
        
    sp = '     '
    for i, line in enumerate(bas[1:]):
        sline = line.split()
        elem = element[sline[0]]
        xyz.append(elem.symbol.rjust(3) + sp + sp.join([f'{float(s):10.5f}' for s in sline[1:]]))

    return read(StringIO('\n'.join(xyz)), format='xyz')

def write_infodat(infodat) -> None:
        string = []
        string.append("   pyreball info.dat")
        string.append(f"  {len(infodat.numbers):2d}\t\t - Number of species")
        for i, num in enumerate(infodat.numbers):
            string.append("  " + 70*"=")
            string.append(f"  {i+1:2d}\t\t - Information for this species")
            string.append(f"  {infodat.elements[num].symbol:<2}\t\t - Element")
            string.append(f"  {num:3d}\t\t - Nuclear Z")
            string.append(
                f"  {infodat.elements[num].atomic_weight:7.3f}\t - Atomic Mass"
            )
            string.append(
                f"  {infodat.numshells[num]:2d}\t\t - Number of shells; "
                "L for each shell"
            )
            string.append(
                "  ".join(["", *[f"{x:1d}" for x in infodat.shells[num]]])
            )
            string.append(
                f"  {infodat.numshellsPP[num]:2d}\t\t - Number of shells; "
                "L for each shell  (Pseudopotential)"
            )
            string.append(
                "  ".join(["", *[f"{x:1d}" for x in infodat.shellsPP[num]]])
            )
            string.append(
                f"  {infodat.cutoffPP[num]:12.4f}\t - Radial cutoffs PP"
            )
            string.append(
                "  ".join(["", *[f"{x:5.2f}" for x in infodat.qneutral[num]]])
            )
            string.append(
                "  ".join(["", *[f"{x:<25}" for x in infodat.wffiles[num]]])
            )
            string.append(
                "  ".join(["", *[f"{x:<25}" for x in infodat.nafiles[num]]])
            )
            string.append(
                f"  {infodat.energy[num]:12.4f}\t - Atomic energy"
            )
            string.append("  " + 70*"=")

        with open("info.dat", "w") as fp:
            fp.write("\n".join(string))
        