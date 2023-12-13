from io import StringIO
from os.path import splitext

from ase.io import read
from ase.io.formats import NameOrFile

from ._zsymbol import z_to_symbol

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

    if len(bas) == (int(bas[0]) + 2):
        xyz = bas
    elif len(bas) == (int(bas[0]) + 1):
        xyz = [bas[0], "no comment\n"]
        for line in bas[1:]:
            xyz.append(line)
    else:
        raise ValueError(f"{filename} is not a valid .bas file")

    for i, line in enumerate(xyz[2:]):
        sline = line.split(' ')
        for j, c in enumerate(sline):
            if c == '':
                sline[j] = ' '
                continue
            try:
                atomz = str(int(c))
                sline[j] = z_to_symbol[atomz]
            except ValueError:
                pass
        xyz[i+2] = ''.join(sline)

    return read(StringIO(''.join(xyz)), format='xyz')