from __future__ import annotations
import pathlib
import numpy as np
import numpy.typing as npt

class WaveFunction:
    def __init__(self, z: int, l: int, q: float, rcut: float, r: npt.NDArray[np.float64], psi: npt.NDArray[np.float64]) -> None:
        self.z = z
        self.l = l
        self.q = q
        self.rcut = rcut
        self.r = r
        self.psi = psi
        if l == 0:
            self.lsym = 's'
        elif l == 1:
            self.lsym = 'p'
        elif l == 2:
            self.lsym = 'd'
        else:
            self.lsym = ''

    @classmethod
    def load(cls, path: pathlib.Path) -> WaveFunction:
        raw = path.read_text().splitlines()
        z, l = list(map(int, raw[0].split('=')[1].split()))
        q = float(raw[1].split('=')[1])
        r, psi = np.loadtxt(raw[3:], unpack=True)
        rcut = r[-1]
        return cls(z, l, q, rcut, r, psi)
