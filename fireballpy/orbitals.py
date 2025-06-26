import numpy as np
from scipy.integrate import tplquad
from ase.data import chemical_symbols
from ase.units import Bohr, Hartree
from numba import njit

from fireballpy.fireball import BaseFireball
from fireballpy.utils import ANGULAR_MOMENTUM

Y00 = 0.5/np.sqrt(np.pi)
Y00SQRT3 = Y00*np.sqrt(3.0)
Y00SQRT5 = Y00*np.sqrt(5.0)


@njit
def lda_exchange_energy(n):
    return -0.7385587664*np.cbrt(n)*Bohr*Hartree

@njit
def lda_exchange_potential(n):
    return -0.9847450219*np.cbrt(n)*Bohr*Hartree

@njit
def lda_exchange_correlation_energy(n):
    ncbrt = np.cbrt(n)*Bohr
    if ncbrt < 1e-16:
        return 0.0
    rs = 0.62035049/ncbrt
    if rs > 1.0:
        sqrs = np.sqrt(rs)
        den = 1.0 + 1.0529*sqrs + 0.3334*rs
        return (-0.4582/rs - 0.1423/den)*Hartree
    lrs = np.log(rs)
    return (-0.4582/rs - 0.0480 + 0.0311*lrs - 0.0116*rs + 0.002*rs*lrs)*Hartree

@njit
def lda_exchange_correlation_potential(n):
    ncbrt = np.cbrt(n)*Bohr
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

@njit
def lda_exchange_correlation_energy_derivative(n):
    if n < 1e-8:
        return 0.0
    exc = lda_exchange_correlation_energy(n)/Hartree
    vxc = lda_exchange_correlation_potential(n)/Hartree
    return Hartree*(vxc - exc)/n

@njit
def lda_exchange_correlation_potential_derivative(n):
    nb = n*Bohr**3
    ncbrt = np.cbrt(nb)
    if ncbrt < 1e-16:
        return 0.0
    ex = lda_exchange_energy(n)/Hartree
    dexc = lda_exchange_correlation_energy_derivative(n)/(Hartree*Bohr**3)
    rs = 0.62035049/ncbrt
    if rs > 1.0:
        sqrs = np.sqrt(rs)
        den = 1.0 + 1.0529*sqrs + 0.3334*rs
        denp =  1.0529/(2.0*sqrs) + 0.3334
        denpp = -0.5*1.0529/(2.0*rs*sqrs)
        depsc = 0.1423*denp/(den*den)
        ddepsc = -2.0*0.1423*denp*denp/den**3 + 0.1423*denpp/(den*den)
    else:
        lrs = np.log(rs)
        depsc = 0.0311/rs - 0.0116 + 0.002*(lrs + 1.0)
        ddepsc = -0.0311/(rs*rs) + 0.002/rs
    dpotc = (4.0*rs/(9.0*nb*nb))*depsc + (rs*rs/(9.0*nb*nb))*ddepsc
    dpotx = -(2.0/(9.0*nb*nb))*ex
    return Hartree*Bohr**3*(2.0*dexc + nb*(dpotx + dpotc))

@njit
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

@njit
def orbital(x, y, z, fcoefs, fxs, rn, l, m):
    dx = x - rn[0]
    dy = y - rn[1]
    dz = z - rn[2]
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    if dr >= np.max(fxs):
        return np.float64(0.0)
    if dr < 1e-10:
        return np.float64(0.0)
    ix = 0
    for i in range(1, fxs.size):
        if dr < fxs[i]:
            ix = i - 1
            break
    dri = dr - fxs[ix]
    fr = fcoefs[3, ix] + dri*(fcoefs[2, ix] + dri*(fcoefs[1, ix] + dri*fcoefs[0, ix]))
    if l == 0:
        return Y00*fr
    fy = sphham(l, m, dx, dy, dz)
    return fr*fy/dr**l

class Orbitals:
    def __init__(self, fbobj: BaseFireball) -> None:
        # Save to self for convenience
        fbobj.fdatafiles.load_wf()
        self.coords = fbobj.atomsystem.positions
        self.wfs = fbobj.fdatafiles.wfs
        self.shell_charges = fbobj.shell_charges
        self.nums = fbobj.atomsystem.numbers

        self.xs = np.zeros((fbobj.natoms, fbobj.maxshells, 1001), dtype=np.float64)
        self.cs = np.zeros((fbobj.natoms, fbobj.maxshells, 4, 1000), dtype=np.float64)
        self.aslm = np.zeros((fbobj.norbitals, 4), dtype=np.int32)
        k = 0
        for i, z in enumerate(self.nums):
            j = 0
            for orbname in ['s', 'p', 'd', 's*', 'p*', 'd*']:
                if orbname not in self.wfs[z]:
                    continue
                self.xs[i, j, :] = self.wfs[z][orbname].x
                self.cs[i, j, :, :] = self.wfs[z][orbname].c
                l = ANGULAR_MOMENTUM[orbname[0]]
                for m in range(-l, l+1):
                    self.aslm[k, 0] = i
                    self.aslm[k, 1] = j
                    self.aslm[k, 2] = l
                    self.aslm[k, 3] = m
                    k += 1
                j += 1

    def get_names(self) -> list[str]:
        names: list[str] = []
        for z in self.nums:
            for orbname in ['s', 'p', 'd', 's*', 'p*', 'd*']:
                if orbname in self.wfs[z]:
                    names.append(chemical_symbols[z] + orbname)
        return names

    def munu_integral(self, mu: int, nu: int, iatom: int = -1, f=None):
        aslm = self.aslm
        scs = self.shell_charges
        rns = self.coords
        fcs = self.cs
        fxs = self.xs
        i1 = aslm[mu][0]
        j1 = aslm[mu][1]
        i2 = aslm[nu][0]
        j2 = aslm[nu][1]
        rn1 = rns[i1]
        rn2 = rns[i2]
        rcut1 = self.xs[i1, j1, :].max()
        rcut2 = self.xs[i2, j2, :].max()
        if i1 == i2:
            x0 = rn1[0]
            y0 = rn1[1]
            z0 = rn1[2]
            xmax = rcut1
            ymax = rcut1
            zmax = rcut1
        else:
            center = 0.5*(rn1 + rn2)
            disth = 0.5*np.abs(rn1 - rn2)
            cutmin = np.minimum(rcut1, rcut2)
            if 2.0*np.min(disth) > cutmin:
                return np.float64(0.0)
            x0 = center[0]
            y0 = center[1]
            z0 = center[2]
            xmax = disth[0] + cutmin
            ymax = disth[1] + cutmin
            zmax = disth[2] + cutmin

        if f is None:
            if iatom > -1:
                @njit
                def ifun0(z, y, x, aiatom, ascs, aaslm, afcs, afxs, arns, mu, nu):
                    i1 = aaslm[mu][0]
                    j1 = aaslm[mu][1]
                    l1 = aaslm[mu][2]
                    m1 = aaslm[mu][3]
                    i2 = aaslm[nu][0]
                    j2 = aaslm[nu][1]
                    l2 = aaslm[nu][2]
                    m2 = aaslm[nu][3]
                    xs1 = afxs[i1, j1, :]
                    xs2 = afxs[i2, j2, :]
                    cs1 = afcs[i1, j1, :, :]
                    cs2 = afcs[i2, j2, :, :]
                    rn1 = arns[i1]
                    rn2 = arns[i2]
                    n = 0.0
                    rn3 = arns[aiatom]
                    for j in range(ascs.shape[1]):
                        q = ascs[aiatom, j]
                        if np.abs(q) < 1e-10:
                            continue
                        xs3 = afxs[aiatom, j, :]
                        cs3 = afcs[aiatom, j, :, :]
                        n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
                    o1 = orbital(x, y, z, cs1, xs1, rn1, l1, m1)
                    o2 = orbital(x, y, z, cs2, xs2, rn2, l2, m2)
                    res = o1*o2*n
                    return res
                return tplquad(ifun0, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                               args=(iatom, scs, aslm, fcs, fxs, rns, mu, nu))[0]

            @njit
            def ifun1(z, y, x, ascs, aaslm, afcs, afxs, arns, mu, nu):
                i1 = aaslm[mu][0]
                j1 = aaslm[mu][1]
                l1 = aaslm[mu][2]
                m1 = aaslm[mu][3]
                i2 = aaslm[nu][0]
                j2 = aaslm[nu][1]
                l2 = aaslm[nu][2]
                m2 = aaslm[nu][3]
                xs1 = afxs[i1, j1, :]
                xs2 = afxs[i2, j2, :]
                cs1 = afcs[i1, j1, :, :]
                cs2 = afcs[i2, j2, :, :]
                rn1 = arns[i1]
                rn2 = arns[i2]
                n = 0.0
                for i in range(ascs.shape[0]):
                    rn3 = arns[i]
                    for j in range(ascs.shape[1]):
                        q = ascs[i, j]
                        if np.abs(q) < 1e-10:
                            continue
                        xs3 = afxs[i, j, :]
                        cs3 = afcs[i, j, :, :]
                        n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
                o1 = orbital(x, y, z, cs1, xs1, rn1, l1, m1)
                o2 = orbital(x, y, z, cs2, xs2, rn2, l2, m2)
                res = o1*o2*n
                return res
            return tplquad(ifun1, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                           args=(scs, aslm, fcs, fxs, rns, mu, nu))[0]

        if isinstance(f, int) or isinstance(f, float):
            @njit
            def ifunn(z, y, x, af, aaslm, afcs, afxs, arns, mu, nu):
                i1 = aaslm[mu][0]
                j1 = aaslm[mu][1]
                l1 = aaslm[mu][2]
                m1 = aaslm[mu][3]
                i2 = aaslm[nu][0]
                j2 = aaslm[nu][1]
                l2 = aaslm[nu][2]
                m2 = aaslm[nu][3]
                xs1 = afxs[i1, j1, :]
                xs2 = afxs[i2, j2, :]
                cs1 = afcs[i1, j1, :, :]
                cs2 = afcs[i2, j2, :, :]
                rn1 = arns[i1]
                rn2 = arns[i2]
                o1 = orbital(x, y, z, cs1, xs1, rn1, l1, m1)
                o2 = orbital(x, y, z, cs2, xs2, rn2, l2, m2)
                res = o1*o2*af
                return res
            return tplquad(ifunn, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                           args=(f, aslm, fcs, fxs, rns, mu, nu))[0]

        if iatom > -1:
            @njit
            def ifun3(z, y, x, aiatom, ascs, aaslm, afcs, afxs, arns, mu, nu):
                i1 = aaslm[mu][0]
                j1 = aaslm[mu][1]
                l1 = aaslm[mu][2]
                m1 = aaslm[mu][3]
                i2 = aaslm[nu][0]
                j2 = aaslm[nu][1]
                l2 = aaslm[nu][2]
                m2 = aaslm[nu][3]
                xs1 = afxs[i1, j1, :]
                xs2 = afxs[i2, j2, :]
                cs1 = afcs[i1, j1, :, :]
                cs2 = afcs[i2, j2, :, :]
                rn1 = arns[i1]
                rn2 = arns[i2]
                n = 0.0
                rn3 = arns[aiatom]
                for j in range(ascs.shape[1]):
                    q = ascs[aiatom, j]
                    if np.abs(q) < 1e-10:
                        continue
                    xs3 = afxs[aiatom, j, :]
                    cs3 = afcs[aiatom, j, :, :]
                    n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
                o1 = orbital(x, y, z, cs1, xs1, rn1, l1, m1)
                o2 = orbital(x, y, z, cs2, xs2, rn2, l2, m2)
                fn = f(n)
                res = o1*o2*fn
                return res
            return tplquad(ifun3, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                           args=(iatom, scs, aslm, fcs, fxs, rns, mu, nu))[0]

        @njit
        def ifun4(z, y, x, ascs, aaslm, afcs, afxs, arns, mu, nu):
            i1 = aaslm[mu][0]
            j1 = aaslm[mu][1]
            l1 = aaslm[mu][2]
            m1 = aaslm[mu][3]
            i2 = aaslm[nu][0]
            j2 = aaslm[nu][1]
            l2 = aaslm[nu][2]
            m2 = aaslm[nu][3]
            xs1 = afxs[i1, j1, :]
            xs2 = afxs[i2, j2, :]
            cs1 = afcs[i1, j1, :, :]
            cs2 = afcs[i2, j2, :, :]
            rn1 = arns[i1]
            rn2 = arns[i2]
            n = 0.0
            for i in range(ascs.shape[0]):
                rn3 = arns[i]
                for j in range(ascs.shape[1]):
                    q = ascs[i, j]
                    if np.abs(q) < 1e-10:
                        continue
                    xs3 = afxs[i, j, :]
                    cs3 = afcs[i, j, :, :]
                    n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
            o1 = orbital(x, y, z, cs1, xs1, rn1, l1, m1)
            o2 = orbital(x, y, z, cs2, xs2, rn2, l2, m2)
            fn = f(n)
            res = o1*o2*fn
            return res
        return tplquad(ifun4, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                       args=(scs, aslm, fcs, fxs, rns, mu, nu))[0]

    def munu_weight_integral(self, mu: int, nu: int, iatom: int = -1, f=None):
        aslm = self.aslm
        scs = self.shell_charges
        rns = self.coords
        fcs = self.cs
        fxs = self.xs
        i1 = aslm[mu][0]
        j1 = aslm[mu][1]
        i2 = aslm[nu][0]
        j2 = aslm[nu][1]
        rn1 = rns[i1]
        rn2 = rns[i2]
        rcut1 = self.xs[i1, j1, :].max()
        rcut2 = self.xs[i2, j2, :].max()
        if i1 == i2:
            x0 = rn1[0]
            y0 = rn1[1]
            z0 = rn1[2]
            xmax = rcut1
            ymax = rcut1
            zmax = rcut1
        else:
            center = 0.5*(rn1 + rn2)
            disth = 0.5*np.abs(rn1 - rn2)
            cutmin = np.minimum(rcut1, rcut2)
            if 2.0*np.min(disth) > cutmin:
                return np.float64(0.0)
            x0 = center[0]
            y0 = center[1]
            z0 = center[2]
            xmax = disth[0] + cutmin
            ymax = disth[1] + cutmin
            zmax = disth[2] + cutmin

        if f is None:
            if iatom > -1:
                @njit
                def ifun0(z, y, x, aiatom, ascs, aaslm, afcs, afxs, arns, mu, nu):
                    i1 = aaslm[mu][0]
                    j1 = aaslm[mu][1]
                    i2 = aaslm[nu][0]
                    j2 = aaslm[nu][1]
                    xs1 = afxs[i1, j1, :]
                    xs2 = afxs[i2, j2, :]
                    cs1 = afcs[i1, j1, :, :]
                    cs2 = afcs[i2, j2, :, :]
                    rn1 = arns[i1]
                    rn2 = arns[i2]
                    n = 0.0
                    rn3 = arns[aiatom]
                    for j in range(ascs.shape[1]):
                        q = ascs[aiatom, j]
                        if np.abs(q) < 1e-10:
                            continue
                        xs3 = afxs[aiatom, j, :]
                        cs3 = afcs[aiatom, j, :, :]
                        n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
                    o1 = orbital(x, y, z, cs1, xs1, rn1, 0, 0)
                    o2 = orbital(x, y, z, cs2, xs2, rn2, 0, 0)
                    res = np.abs(o1*o2)*n
                    return res
                return tplquad(ifun0, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                               args=(iatom, scs, aslm, fcs, fxs, rns, mu, nu))[0]

            @njit
            def ifun1(z, y, x, ascs, aaslm, afcs, afxs, arns, mu, nu):
                i1 = aaslm[mu][0]
                j1 = aaslm[mu][1]
                i2 = aaslm[nu][0]
                j2 = aaslm[nu][1]
                xs1 = afxs[i1, j1, :]
                xs2 = afxs[i2, j2, :]
                cs1 = afcs[i1, j1, :, :]
                cs2 = afcs[i2, j2, :, :]
                rn1 = arns[i1]
                rn2 = arns[i2]
                n = 0.0
                for i in range(ascs.shape[0]):
                    rn3 = arns[i]
                    for j in range(ascs.shape[1]):
                        q = ascs[i, j]
                        if np.abs(q) < 1e-10:
                            continue
                        xs3 = afxs[i, j, :]
                        cs3 = afcs[i, j, :, :]
                        n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
                o1 = orbital(x, y, z, cs1, xs1, rn1, 0, 0)
                o2 = orbital(x, y, z, cs2, xs2, rn2, 0, 0)
                res = np.abs(o1*o2)*n
                return res
            return tplquad(ifun1, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                           args=(scs, aslm, fcs, fxs, rns, mu, nu))[0]

        if isinstance(f, int) or isinstance(f, float):
            @njit
            def ifunn(z, y, x, af, aaslm, afcs, afxs, arns, mu, nu):
                i1 = aaslm[mu][0]
                j1 = aaslm[mu][1]
                i2 = aaslm[nu][0]
                j2 = aaslm[nu][1]
                xs1 = afxs[i1, j1, :]
                xs2 = afxs[i2, j2, :]
                cs1 = afcs[i1, j1, :, :]
                cs2 = afcs[i2, j2, :, :]
                rn1 = arns[i1]
                rn2 = arns[i2]
                o1 = orbital(x, y, z, cs1, xs1, rn1, 0, 0)
                o2 = orbital(x, y, z, cs2, xs2, rn2, 0, 0)
                res = np.abs(o1*o2)*af
                return res
            return tplquad(ifunn, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                           args=(f, aslm, fcs, fxs, rns, mu, nu))[0]

        if iatom > -1:
            @njit
            def ifun3(z, y, x, aiatom, ascs, aaslm, afcs, afxs, arns, mu, nu):
                i1 = aaslm[mu][0]
                j1 = aaslm[mu][1]
                i2 = aaslm[nu][0]
                j2 = aaslm[nu][1]
                xs1 = afxs[i1, j1, :]
                xs2 = afxs[i2, j2, :]
                cs1 = afcs[i1, j1, :, :]
                cs2 = afcs[i2, j2, :, :]
                rn1 = arns[i1]
                rn2 = arns[i2]
                n = 0.0
                rn3 = arns[aiatom]
                for j in range(ascs.shape[1]):
                    q = ascs[aiatom, j]
                    if np.abs(q) < 1e-10:
                        continue
                    xs3 = afxs[aiatom, j, :]
                    cs3 = afcs[aiatom, j, :, :]
                    n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
                o1 = orbital(x, y, z, cs1, xs1, rn1, 0, 0)
                o2 = orbital(x, y, z, cs2, xs2, rn2, 0, 0)
                fn = f(n)
                res = np.abs(o1*o2)*fn
                return res
            return tplquad(ifun3, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                           args=(iatom, scs, aslm, fcs, fxs, rns, mu, nu))[0]

        @njit
        def ifun4(z, y, x, ascs, aaslm, afcs, afxs, arns, mu, nu):
            i1 = aaslm[mu][0]
            j1 = aaslm[mu][1]
            i2 = aaslm[nu][0]
            j2 = aaslm[nu][1]
            xs1 = afxs[i1, j1, :]
            xs2 = afxs[i2, j2, :]
            cs1 = afcs[i1, j1, :, :]
            cs2 = afcs[i2, j2, :, :]
            rn1 = arns[i1]
            rn2 = arns[i2]
            n = 0.0
            for i in range(ascs.shape[0]):
                rn3 = arns[i]
                for j in range(ascs.shape[1]):
                    q = ascs[i, j]
                    if np.abs(q) < 1e-10:
                        continue
                    xs3 = afxs[i, j, :]
                    cs3 = afcs[i, j, :, :]
                    n += q*orbital(x, y, z, cs3, xs3, rn3, 0, 0)**2
            o1 = orbital(x, y, z, cs1, xs1, rn1, 0, 0)
            o2 = orbital(x, y, z, cs2, xs2, rn2, 0, 0)
            fn = f(n)
            res = np.abs(o1*o2)*fn
            return res
        return tplquad(ifun4, x0-xmax, x0+xmax, y0-ymax, y0+ymax, z0-zmax, z0+zmax, epsabs=1e-4, epsrel=1e-4,
                       args=(scs, aslm, fcs, fxs, rns, mu, nu))[0]
