from dataclasses import dataclass, field
from itertools import product
from os.path import join
import numpy as np

from infodat import InfoDat
from utils import (
    alloc_float,
    alloc_int,
    load_data,
    read_line,
    read_float_array,
    read_int_entry,
    skip_lines
)


_NFOFX = 207  # AQUI: que demonios es esto???

_TWOCFNAMES = [
    'overlap',
    'vna_ontopl',
    'vna_ontopr',
    'vna_atom',
    'vnl',
    'xc_ontop',
    'xc_atom',
    'xc_corr',
    'dipole_z',
    'dipole_y',
    'dipole_x',
    'coulomb',
    'kinetic',
    'nuxc',
    'den_ontopl',
    'den_ontopr',
    'den_atom',
    'dnuxc_ol',
    'dnuxc_or',
    'denS_ontopl',
    'denS_ontopr',
    'denS_atom',
    'overlapS',
]


# Useful wrapper
def _penum(it: list, rep: int = 1):
    return product(enumerate(it, repeat=rep))


@dataclass
class _MuNu:
    mu: np.ndarray[int, int]
    nu: np.ndarray[int, int]
    mval: np.ndarray[int, int]


@dataclass
class _OneCenter:
    exc1: np.ndarray[int, float] = field(init=False)
    nuxc: np.ndarray[int, float] = field(init=False)
    dexc: np.ndarray[int, float] = field(init=False)
    dnuxc: np.ndarray[int, float] = field(init=False)


@dataclass
class _TwoCenter:
    # Sizes
    me_max: int = field(init=False)
    me_pp_max: int = field(init=False)
    iter_max: int = field(init=False)

    # Indexes
    index_max: np.ndarray[int, int] = field(init=False)

    # Data
    numz: np.ndarray[int, int] = field(init=False)
    cl_pp: np.ndarray[int, float] = field(init=False)
    zmax: np.ndarray[int, float] = field(init=False)
    xint: np.ndarray[int, float] = field(init=False)


@dataclass
class _ThreeCenter:
    # Sizes
    me_max: int = field(init=False)
    me_pp_max: int = field(init=False)

    # Data
    nbcna: np.ndarray[int, int] = field(init=False)
    mbcna: np.ndarray[int, float] = field(init=False)
    hbcna: np.ndarray[int, float] = field(init=False)
    bcna: np.ndarray[int, float] = field(init=False)
    nden: np.ndarray[int, int] = field(init=False)
    mden: np.ndarray[int, float] = field(init=False)
    hden: np.ndarray[int, float] = field(init=False)
    den: np.ndarray[int, float] = field(init=False)
    dens: np.ndarray[int, float] = field(init=False)


@dataclass
class _FData:
    nums: list[int]
    fpath: str
    idat: InfoDat

    # Derived quantities
    nspecies: int = field(init=False)
    nsh_max: int = field(init=False)
    norb: np.ndarray[int, int] = field(init=False)
    norb_pp: np.ndarray[int, int] = field(init=False)

    # Subclasses for arrays
    munu: _MuNu = field(default_factory=_MuNu, repr=False)
    onec: _OneCenter = field(default_factory=_OneCenter, repr=False)
    twoc: _TwoCenter = field(default_factory=_TwoCenter, repr=False)
    threec: _ThreeCenter = field(default_factory=_ThreeCenter, repr=False)

    def __post_init__(self):
        # Derive quantities
        self.nspecies = len(self.nums)
        self.nsh_max = max(
            max(len(l) for _, l in self.idat.shells.items()),
            max(len(l) for _, l in self.idat.shells_pp.items()))
        self.norbs = np.array(
            [np.sum(2*self.idat.shells[num] + 1) for num in self.nums],
            dtype=int, order="F")

        # Create auxiliary arrays
        self._get_maxdims()
        self._make_munu()

        # Load all data
        self._load_one_center()
        self._load_two_center()
        self._load_three_center()

    def _get_maxdims(self) -> None:
        me2c_max, me3c_max = 0, 0
        for (in1, num1), (in2, num2) in _penum(self.nums, 2):
            idx = 0
            for l1, l2 in product(self.idat.shells[num1],
                                  self.idat.shells[num2]):
                idx += 2*min(l1, l2) + 1
            me2c_max = max(idx, me2c_max)
            for l1, l2 in product(self.idat.shells[num1],
                                  self.idat.shells[num2]):
                if l1 == l2 == 0:
                    continue
                l2z, l22 = int(l2 != 0), int(l2 == 2)
                idx += l22 + (l1 == 0) + \
                    (l1 == 1)*(1 + l2z + 2*l22) + \
                    (l1 == 2)*(1 + 3*l2z + 2*l22)
            me3c_max = max(idx, me3c_max)

        self.twoc.me_max = me2c_max
        self.threec.me_pp_max = me3c_max

    def _make_munu(self) -> None:
        mu = alloc_float(self.threec.me_max, self.nspecies, self.nspecies)
        nu = alloc_float(self.threec.me_max, self.nspecies, self.nspecies)
        mval = alloc_float(self.threec.me_max, self.nspecies, self.nspecies)
        for (in1, num1), (in2, num2) in _penum(self.nums, 2):
            idx = 0
            n1 = 0
            for l1 in self.idat.shells[num1]:
                n1 += l1 + 1
                n2 = 0
                for l2 in self.idat.shells[num2]:
                    n2 += l2 + 1
                    for imu in range(min(l1, l2), min(l1, l2) + 1):
                        mu[idx, in1, in2] = n1 + imu
                        nu[idx, in1, in2] = n2 + imu
                        idx += 1
                    n2 += l2
                n1 += l1
            self.twoc.im[in1, in2] = idx

            # AQUI: l > 2 posible ??
            n1 = 0
            for l1 in self.idat.shells[num1]:
                n1 += l1 + 1
                n2 = 0
                for l2 in self.idat.shells[num2]:
                    n2 += l2 + 1
                    if l1 == l2 == 0:
                        continue
                    mu[idx, in1, in2] = n1 + (l1 > 0)
                    nu[idx, in1, in2] = n2 + (l1 == 0)
                    mval[idx, in1, in2] = 1
                    idx += 1
                    if l2 == 0:
                        continue
                    mu[idx, in1, in2] = n1
                    nu[idx, in1, in2] = n2 + 1
                    mval[idx, in1, in2] = 1
                    idx += 1
                    if l1 == 2:
                        mu[idx, in1, in2] = n1 - 2
                        nu[idx, in1, in2] = n2 - 1
                        mval[idx, in1, in2] = 1
                        idx += 1
                        mu[idx, in1, in2] = n1 + 2
                        nu[idx, in1, in2] = n2 + 1
                        mval[idx, in1, in2] = 1
                        idx += 1
                    if l2 == 2:
                        mu[idx, in1, in2] = n1 + 1
                        nu[idx, in1, in2] = n2 + 2
                        mval[idx, in1, in2] = 1
                        idx += 1
                        mu[idx, in1, in2] = n1 - 1
                        nu[idx, in1, in2] = n2 - 2
                        mval[idx, in1, in2] = 1
                        idx += 1
                    n2 += l2
                n1 += l1

            n1 = 0
            for l1 in self.idat.shells[num1]:
                n1 += l1 + 1
                n2 = 0
                for l2 in self.idat.shells[num2]:
                    n2 += l2 + 1
                    if l1 == 2:
                        mu[idx, in1, in2] = n1 + 2
                        nu[idx, in1, in2] = n2
                        mval[idx, in1, in2] = 2
                        idx += 1
                    if l2 == 2:
                        mu[idx, in1, in2] = n1
                        nu[idx, in1, in2] = n2 + 2
                        mval[idx, in1, in2] = 2
                        idx += 1
                    n2 += l2
                n1 += l1
            self.threec.im[in1, in2] = idx

    def _load_one_center(self) -> None:
        exc = alloc_float(self.nspecies, self.nsh_max, self.nsh_max)
        nuxc = alloc_float(self.nspecies, self.nsh_max, self.nsh_max)
        dexc = alloc_float(self.nspecies, self.nsh_max,
                           self.nsh_max, self.nsh_max)
        dnuxc = alloc_float(self.nspecies, self.nsh_max,
                            self.nsh_max, self.nsh_max)

        for i, num in enumerate(self.nums):
            # exc and nuxc
            data = load_data(
                join(self.fpath, f"xc1c_dq1.{num:02}.dat"),
                header=6)
            niter = read_int_entry(data, field=1)
            exc[i, :niter, :niter] = read_float_array(
                data, lines=niter, width=niter)
            skip_lines(data)
            nuxc[i, :niter, :niter] = read_float_array(
                data, lines=niter, width=niter)

            # dexc
            data = load_data(
                join(self.fpath, f"exc1crho.{num:02}.dat"),
                header=6)
            for j in range(self.idat.numshells[num]):
                niter = read_int_entry(data, field=1)
                dexc[i, :niter, :niter, j] = read_float_array(
                    data, lines=niter, width=niter)

            # dnuxc
            data = load_data(
                join(self.fpath, f"nuxc1crho.{num:02}.dat"),
                header=6)
            for j in range(self.idat.numshells[num]):
                niter = read_int_entry(data, field=1)
                dnuxc[i, :niter, :niter, j] = read_float_array(
                    data, lines=niter, width=niter)

        self.onec.exc = exc
        self.onec.nuxc = nuxc
        self.onec.dexc = dexc
        self.onec.dnuxc = dnuxc

    def _load_two_center(self) -> None:
        numz = alloc_int(self.twoc.iter_max, self.nspecies, self.nspecies)
        zmax = alloc_float(self.twoc.iter_max, self.nspecies, self.nspecies)
        xint = alloc_float(self.me_max, _NFOFX, self.twoc.iter_max,
                           self.nspecies, self.nspecies)
        cl_pp = alloc_float(self.nsh_max, self.nspecies)

        for (in1, num1), (in2, num2) in _penum(self.nums, 2):
            for (i, fname) in enumerate(_TWOCFNAMES):
                if i == 13:  # nuxc still to be added
                    continue

                match i:
                    case 1 | 14 | 17 | 19:
                        isub = self.idat.numshells[num1]
                    case 2 | 3 | 15 | 16 | 18 | 20 | 21:
                        isub = self.idat.numshells[num2]
                    case 5 | 6 | 7:
                        isub = 4
                    case _:
                        isub = 0

                for isorp in range(isub + 1):
                    if (i > 13) and (i != 22):
                        continue
                    ffname = f"{fname}_{isorp:02}" if isorp else fname
                    ffname += f".{num1:02}.{num2:02}.dat"
                    data = load_data(join(self.fpath, ffname), header=11)
                    if i == 4:
                        npp = read_int_entry(data)
                        cl_pp[:npp, in1] = read_float_array(data, width=npp)
                    zm, nz = read_line(data, "float", "int")
                    if nz > _NFOFX:
                        raise EnvironmentError("numz > nfofx")
                    zmax[self.twoc.index[i, isorp], in1, in2] = zm
                    numz[self.towc.index[i, isorp], in1, in2] = nz

                    match i:
                        case 3 | 6 | 16:
                            nnz = self.twoc.index_max[in1, in1]
                        case 4:
                            nnz = self.twoc.index_pp_max[in1, in2]
                        case 7:
                            nnz = 1
                        case 9:
                            nnz = self.twoc.index_dy_max[in1, in2]
                        case 10:
                            nnz = self.twoc.index_dx_max[in1, in2]
                        case 11 | 13:
                            nnz = self.idat.numshells[num1] * \
                                self.idat.numshells[num2]
                        case 19 | 20 | 22:
                            nnz = self.index_s_max[in1, in2]
                        case 21:
                            nnz = self.index_s_max[in1, in1]
                        case _:
                            nnz = self.twoc.index_max[in1, in2]

                    xint[:nnz, :nz, self.twoc.index[i, isorp], in1, in2] = \
                        read_float_array(data, lines=nz, width=nnz) \
                        .transpose()

        self.twoc.numz = numz
        self.twoc.cl_pp = cl_pp
        self.twoc.zmax = zmax
        self.twoc.xint = xint

    def _load_three_center(self) -> None:
        self.threec.map = np.array(
            [[[k + (j + i*self.nspecies)*self.nspecies + 1
             for k in range(self.nspecies)]
             for j in range(self.nspecies)]
             for i in range(self.nspecies)],
            dtype=int, order="F")

        # Interaction 1
        nis = alloc_int(2, self.nsh_max + 1, self.nspecies**3)
        mis = alloc_float(2, self.nsh_max + 1, self.nspecies**3)
        his = alloc_float(2, self.nsh_max + 1, self.nspecies**3)
        ais = alloc_float(self.threec.nxm, self.threec.nym, self.threec.me_max,
                          self.nsh_max + 1, self.nspecies**3, 5)

        for it in range(5):
            for (in1, num1), (in2, num2), (in3, num3) in _penum(self.nums, 3):
                idx = self.threec.map[in1, in2, in3] - 1
                nnz = self.threec.index_max[in1, in2]
                for isorp in range(self.idat.numshells[num3] + 1):
                    ffname = f"bcna_{it:02}_{isorp:02}"
                    ffname += f".{num1:02}.{num2:02}.{num3:02}.dat"
                    data = load_data(join(self.fpath, ffname), header=11)
                    ym, ny = read_line(data, "float", "int")
                    xm, nx = read_line(data, "float", "int")
                    skip_lines(data, lines=5)

                    nis[:, isorp, idx] = (nx, ny)
                    mis[:, isorp, idx] = (xm, ym)
                    his[0, isorp, idx] = xm / (nx - 1)
                    his[1, isorp, idx] = ym / (ny - 1)
                    for iy in range(ny):
                        ais[:nx, ny, :nnz, isorp, idx, it] = \
                            read_float_array(data, lines=nx, width=nnz) \
                            .transpose()

        self.threec.nbcna = nis
        self.threec.mbcna = mis
        self.threec.hbcna = his
        self.threec.bcna = ais

        # Interaction 3
        nis = alloc_int(2, self.nsh_max + 1, self.nspecies**3)
        mis = alloc_float(2, self.nsh_max + 1, self.nspecies**3)
        his = alloc_float(2, self.nsh_max + 1, self.nspecies**3)
        ais = alloc_float(self.threec.nxm, self.threec.nym, self.threec.me_max,
                          self.nsh_max + 1, self.nspecies**3, 5)

        for it in range(5):
            for (in3, num3), (in2, num2), (in3, num3) in _penum(self.nums, 3):
                idx = self.threec.map[in1, in2, in3] - 1
                nnz = self.threec.index_max[in1, in2]
                for isorp, it in range(1, self.idat.numshells[num3] + 1):
                    ffname = f"den3_{it:02}_{isorp:02}"
                    ffname += f".{num1:02}.{num2:02}.{num3:02}.dat"
                    data = load_data(join(self.fpath, ffname), header=11)
                    ym, ny = read_line(data, "float", "int")
                    xm, nx = read_line(data, "float", "int")
                    skip_lines(data, lines=5)

                    nis[:, isorp, idx] = (nx, ny)
                    mis[:, isorp, idx] = (xm, ym)
                    his[0, isorp, idx] = xm / (nx - 1)
                    his[1, isorp, idx] = ym / (ny - 1)
                    for iy in range(ny):
                        ais[:nx, ny, :nnz, isorp, idx, it] = \
                            read_float_array(data, lines=nx, width=nnz) \
                            .transpose()

        self.threec.nden = nis
        self.threec.mden = mis
        self.threec.hden = his
        self.threec.den = ais

        # Interaction 4
        ais = alloc_float(self.threec.nxm, self.threec.nym, self.threec.me_max,
                          self.nsh_max + 1, self.nspecies**3, 5)

        for it in range(5):
            for (in1, num1), (in2, num2), (in3, num3) in _penum(self.nums, 3):
                idx = self.threec.icon[in1, in2, in3] - 1
                nnz = self.index_s_max[in1, in2]
                for isorp in range(1, self.idat.numshells[num3] + 1):
                    ffname = f"deS3_{it:02}_{isorp:02}"
                    ffname += f".{num1:02}.{num2:02}.{num3:02}.dat"
                    data = load_data(join(self.fpath, ffname), header=11)
                    ym, ny = read_line(data, "float", "int")
                    xm, nx = read_line(data, "float", "int")
                    skip_lines(data, lines=5)

                    for iy in range(ny):
                        ais[:nx, ny, :nnz, isorp, idx, it] = \
                            read_float_array(data, lines=nx, width=nnz) \
                            .transpose()

        self.threec.dens = ais
