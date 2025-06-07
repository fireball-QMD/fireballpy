import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import PPoly
from scipy.integrate import tplquad
from ase.data import chemical_symbols
from ase.units import Bohr, Hartree

from fireballpy.atoms import AtomSystem
from fireballpy.fdata import FDataFiles
from fireballpy.utils import ANGULAR_MOMENTUM

Y00 = 0.5/np.sqrt(np.pi)
Y00SQRT3 = Y00*np.sqrt(3.0)
Y00SQRT5 = Y00*np.sqrt(5.0)


def lda_exchange_energy(n):
    def ex(x, y, z):
        return -0.7385587664*np.cbrt(n(x,y,z))*Bohr*Hartree
    return ex

def lda_exchange_potential(n):
    def fx(x, y, z):
        return -0.9847450219*np.cbrt(n(x,y,z))*Bohr*Hartree
    return fx

def lda_exchange_correlation_energy(n):
    def exc(x, y, z):
        ncbrt = np.cbrt(n(x, y, z))*Bohr
        if ncbrt < 1e-16:
            return 0.0
        rs = 0.62035049/ncbrt
        if rs > 1.0:
            sqrs = np.sqrt(rs)
            den = 1.0 + 1.0529*sqrs + 0.3334*rs
            return (-0.4582/rs - 0.1423/den)*Hartree
        lrs = np.log(rs)
        return (-0.4582/rs - 0.0480 + 0.0311*lrs - 0.0116*rs + 0.002*rs*lrs)*Hartree
    return exc

def lda_exchange_correlation_potential(n):
    def fxc(x, y, z):
        ncbrt = np.cbrt(n(x, y, z))*Bohr
        if ncbrt < 1e-16:
            return 0.0
        rs = 0.62035049/ncbrt
        if rs > 1.0:
            sqrs = np.sqrt(rs)
            den = 1.0 + 1.0529*sqrs + 0.3334*rs
            return (-0.4582/rs - 0.1423/den - rs*(0.15273333/rs**2 + (0.02497128/sqrs + 0.01581427)/den**2))*Hartree
        lrs = np.log(rs)
        return (-0.4582/rs - 0.0480 + 0.0311*lrs - 0.0116*rs + 0.002*rs*lrs \
            - rs*(0.15273333/rs**2 + 0.01036667/rs - 0.003866667 + 0.00066667*(1.0 + lrs)))*Hartree
    return fxc


def sphham(l: int, m: int, x: np.float64, y: np.float64, z: np.float64) -> np.float64:
    """Factor r^(-l) not included"""
    if l == 0:
      return Y00
    if l == 1:
      if m == -1: sph = y
      elif m ==  1: sph = x
      else: sph = z
      return sph*Y00SQRT3
    if m == -2: sph = 1.7320508075688772*x*y
    elif m ==  2: sph = 0.8660254037844386*(x*x - y*y)
    elif m == -1: sph = 1.7320508075688772*y*z
    elif m ==  1: sph = 1.7320508075688772*x*z
    else: sph = 0.5*(2*z*z - x*x - y*y)
    return sph*Y00SQRT5


class Orbital:
    def __call__(self, *args) -> np.float64:
        if len(args) == 1:
            r = args[0]
        elif len(args) == 3:
            r = np.array(args, dtype=np.float64)
        else:
            raise ValueError("Incorrect number of inputs")
        d = r - self.rn
        dr = np.sqrt(np.sum(d**2))
        if dr >= self.rcutoff:
            return np.float64(0.0)
        fr = np.float64(self.f(dr))
        if self.m is None or self.l == 0:
            return Y00*fr
        if dr < 1e-10:
            return np.float64(0.0)
        fy = sphham(self.l, self.m, d[0], d[1], d[2])
        return fr*fy/dr**self.l

    def __init__(self, f: PPoly, iatom:int, rn: NDArray[np.float64], l: int, m: int | None = None):
        self.f = f
        self.iatom = iatom
        self.l = l
        self.m = m
        self.rn = rn
        self.rcutoff = np.max(f.x)


class Orbitals:
    def __init__(self, atoms: AtomSystem, fdata: FDataFiles) -> None:
        # Save to self for convenience
        self.atoms = atoms
        self.fdata = fdata
        self.fdata.load_wf()
        self.wfs = self.fdata.wfs

        # Two useful indexes
        self.orbitals_map: list[list[int]] = []
        self.shells_map: list[list[int]] = []
        self.orbitals: list[Orbital] = []
        self.shells: list[Orbital] = []
        j = 0
        k = 0
        for i, (z, p) in enumerate(zip(self.atoms.numbers, self.atoms.positions)):
            self.orbitals_map.append([])
            self.shells_map.append([])
            for orbname in ['s', 'p', 'd', 's*', 'p*', 'd*']:
                if orbname not in self.wfs[z]:
                    continue
                l = ANGULAR_MOMENTUM[orbname[0]]
                self.shells.append(Orbital(self.wfs[z][orbname], i, p, l))
                self.shells_map[-1].append(j)
                j += 1
                for m in range(-l, l+1):
                    self.orbitals.append(Orbital(self.wfs[z][orbname], i, p, l, m))
                    self.orbitals_map[-1].append(k)
                    k += 1

    def get_names(self) -> list[str]:
        names: list[str] = []
        for z in self.atoms.numbers:
            for orbname in ['s', 'p', 'd', 's*', 'p*', 'd*']:
                if orbname in self.wfs[z]:
                    names.append(chemical_symbols[z] + orbname)
        return names

    def sph_density(self, shell_charges: NDArray[np.float64]):
        def n(x, y, z):
            res = 0.0
            for i, qs in enumerate(shell_charges):
                for j, q in enumerate(qs):
                    if np.abs(q) < 1e-8:
                        continue
                    res += q*self.shells[self.shells_map[i][j]](x, y, z)**2
            return res
        return n

    def sph_density_atom(self, iatom:int, shell_charges: NDArray[np.float64]):
        def ni(x, y, z):
            res = 0.0
            for j, q in enumerate(shell_charges[iatom]):
                if np.abs(q) < 1e-8:
                    continue
                res += q*self.shells[self.shells_map[iatom][j]](x, y, z)**2
            return res
        return ni

    def munu_integral(self, f, mu: int, nu: int):
        phimu = self.orbitals[mu]
        phinu = self.orbitals[nu]
        if phimu.iatom == phinu.iatom:
            xmin = phimu.rn[0] - phimu.rcutoff
            ymin = phimu.rn[1] - phimu.rcutoff
            zmin = phimu.rn[2] - phimu.rcutoff
            xmax = phimu.rn[0] + phimu.rcutoff
            ymax = phimu.rn[1] + phimu.rcutoff
            zmax = phimu.rn[2] + phimu.rcutoff
            return tplquad(lambda z, y, x: phimu(x, y, z)*f(x, y, z)*phinu(x, y, z),
                           xmin, xmax, ymin, ymax, zmin, zmax, epsabs=1e-2, epsrel=1e-2)
        center = 0.5*(phimu.rn + phinu.rn)
        disth = 0.5*np.abs(phimu.rn - phinu.rn)
        cutmin = np.minimum(phimu.rcutoff, phinu.rcutoff)
        if 2.0*np.min(disth) > cutmin:
            return np.float64(0.0), np.float64(0.0)
        xmin = center[0] - disth[0] - cutmin
        ymin = center[1] - disth[1] - cutmin
        zmin = center[2] - disth[2] - cutmin
        xmax = center[0] + disth[0] + cutmin
        ymax = center[1] + disth[1] + cutmin
        zmax = center[2] + disth[2] + cutmin
        return tplquad(lambda z, y, x: phimu(x, y, z)*f(x, y, z)*phinu(x, y, z),
                       xmin, xmax, ymin, ymax, zmin, zmax, epsabs=1e-2, epsrel=1e-2)
