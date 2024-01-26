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