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


@dataclass
class _OneCenter:
    exc1: np.ndarray[int, float] = field(init=False)
    nuxc: np.ndarray[int, float] = field(init=False)
    dexc: np.ndarray[int, float] = field(init=False)
    dnuxc: np.ndarray[int, float] = field(init=False)


@dataclass
class _TwoCenter:
    numz: np.ndarray[int, int] = field(init=False)
    cl_pp: np.ndarray[int, float] = field(init=False)
    zmax: np.ndarray[int, float] = field(init=False)
    xint: np.ndarray[int, float] = field(init=False)


@dataclass
class _ThreeCenter:
    nbcna: np.ndarray[int, int] = field(init=False)
    mbcna: np.ndarray[int, float] = field(init=False)
    bcna: np.ndarray[int, float] = field(init=False)
    nden: np.ndarray[int, int] = field(init=False)
    mden: np.ndarray[int, float] = field(init=False)
    den: np.ndarray[int, float] = field(init=False)
    dens: np.ndarray[int, float] = field(init=False)


@dataclass
class _FData:
    nsh_max: int
    nshPP_max: int
    isorpmax: int
    isorpmax_xc: int
    nspecies: int

    one_center: _OneCenter = field(init=False)
    two_center: _TwoCenter = field(init=False)
    three_center: _ThreeCenter = field(init=False)

    def load_1c(self, infodat: InfoDat, fpath: str) -> None:
        # Init arrays
        exc = alloc_float(self.nsh_max, self.nsh_max, self.nspecies)
        nuxc = alloc_float(self.nsh_max, self.nsh_max, self.nspecies)
        dexc = alloc_float(self.nsh_max, self.nsh_max,
                           self.nsh_max, self.nspecies)
        dnuxc = alloc_float(self.nsh_max, self.nsh_max,
                            self.nsh_max, self.nspecies)

        for i, num in enumerate(infodat.numbers):
            nsh = infodat.numshells[num]

            # exc & nuxc
            data = load_data(join(fpath, f"xc1c_dq1.{num:02}.dat"), header=6)
            niter = read_int_entry(data, field=1)
            exc[:niter, :niter, i] = read_float_array(data, lines=niter)
            skip_lines(data)
            nuxc[:niter, :niter, i] = read_float_array(data, lines=niter)

            # dexc
            data = load_data(join(fpath, f"exc1crho.{num:02}.dat"), header=6)
            for j in range(nsh):
                niter = read_int_entry(data, field=1)
                dexc[:niter, :niter, j, i] = read_float_array(
                    data, lines=niter)

            # dnuxc
            data = load_data(join(fpath, f"nuxc1crho.{num:02}.dat"), header=6)
            for j in range(nsh):
                niter = read_int_entry(data, field=1)
                dnuxc[:niter, :niter, j, i] = read_float_array(
                    data, lines=niter)

        self.one_center.exc = exc
        self.one_center.nuxc = nuxc
        self.one_center.dexc = dexc
        self.one_center.dnuxc = dnuxc

        def _two_center(self, infodat: InfoDat, fpath: str) -> None:
            numz = alloc_int(self.i2c_max, self.nspecies, self.nspecies)
            cl_pp = alloc_float(self.nshPP_max, self.nspecies)
            zmax = alloc_float(self.nshPP_max, self.nspecies)
            xint = alloc_float(self.me2c_max, _NFOFX, self.i2c_max,
                               self.nspecies, self.nspecies)

            for i, fname in enumerate(_TWOCFNAMES):
                if i == 13:  # nuxc still to be added
                    continue

                for (in1, num1), (in2, num2) in product(
                        enumerate(infodat.numbers),
                        enumerate(infodat.numbers)):
                    match i:
                        case 1 | 14 | 17 | 19:
                            isub2c = infodat.numshells[num1]
                        case 2 | 3 | 15 | 16 | 18 | 20 | 21:
                            isub2c = infodat.numshells[num2]
                        case 5 | 6 | 7:
                            isub2c = 4
                        case _:
                            isub2c = 0

                    for isorp in range(isub2c + 1):
                        if (i > 13) and (i != 22):
                            continue
                        ffname = f"{fname}_{isorp:02}" if isorp else fname
                        ffname += f".{num1:02}.{num2:02}.dat"
                        data = load_data(join(fpath, ffname), header=11)
                        if i == 4:
                            npp = read_int_entry(data)
                            cl_pp[:npp, in2] = read_float_array(data)
                        zm, nz = read_line(data, "float", "int")
                        if nz > _NFOFX:
                            raise EnvironmentError("numz > nfofx")
                        zmax[self.ind2c[i, isorp], in1, in2] = zm
                        numz[self.ind2c[i, isorp], in1, in2] = nz

                        match i:
                            case 3 | 6:
                                nnz = self.ind2c[in1, in1]
                            case 4:
                                nnz = self.indmaxPP[in1, in2]
                            case 7:
                                nnz = 1
                            case 9:
                                nnz = self.indmax2cDY[in1, in2]
                            case 10:
                                nnz = self.indmax2cDX[in1, in2]
                            case 11 | 13:
                                nnz = num1 * num2
                            case 16:
                                nnz = self.indmax2c[in1, in1]
                            case 19 | 20 | 22:
                                nnz = self.indmaxS[in1, in2]
                            case 21:
                                nnz = self.indmaxS[in1, in1]
                            case _:
                                nnz = self.ind2c[i, isorp]

                        xint[:nnz, :nz, self.ind2c[i, isorp], in1, in2] = \
                            read_float_array(data, lines=nz) \
                            .transpose()[:nnz, :]

            self.two_center.numz = numz
            self.two_center.cl_pp = cl_pp
            self.two_center.zmax = zmax
            self.two_center.xint = xint

    def _three_center(self, infodat: InfoDat, fpath: str) -> None:
        self.three_center.icon = np.array(
            [[[k + (j + i*self.nspecies)*self.nspecies + 1
             for k in range(self.nspecies)]
             for j in range(self.nspecies)]
             for i in range(self.nspecies)],
            dtype=int, order="F")

        # Interactions 1
        nbcna = alloc_int(self.isorpmax + 1, self.nspecies**3, 2)
        mbcna = alloc_float(self.isorpmax + 1, self.nspecies**3, 2)
        bcna = alloc_float(self.me3c_max, self.nxm, self.nym, 5,
                           self.isorpmax + 1, self.nspecies**3)

        for (in3, num3), (in2, num2), (in1, num1) in product(
                enumerate(infodat.numbers),
                enumerate(infodat.numbers),
                enumerate(infodat.numbers)):
            idx = self.three_center.icon[in1, in2, in3] - 1
            nnz = self.three_center.im[in1, in2]
            for isorp, it in product(
                    range(infodat.numshells[num3] + 1),
                    range(5)):
                ffname = f"bcna_{it:02}_{isorp:02}"
                ffname += f".{num1:02}.{num2:02}.{num3:02}.dat"
                data = load_data(join(fpath, ffname), header=11)
                ym, ny = read_line(data, "float", "int")
                xm, nx = read_line(data, "float", "int")
                skip_lines(data, lines=5)

                mbcna[1, isorp, idx], nbcna[1, isorp, idx] = ym, ny
                mbcna[0, isorp, idx], nbcna[0, isorp, idx] = xm, nx
                for iy in range(ny):
                    bcna[:nnz, :nx, iy, it, isorp, idx] = \
                        read_float_array(data, lines=nx).transpose()[:nnz, :]

        self.three_center.nbcna = nbcna
        self.three_center.mbcna = mbcna
        self.three_center.bcna = bcna

        # Interaction 3
        nden = alloc_int(self.isorpmax_xc + 1, self.nspecies**3, 2)
        mden = alloc_float(self.isorpmax_xc + 1, self.nspecies**3, 2)
        den = alloc_float(self.nxm, self.nym, self.me3c_max,
                          self.isorpmax_xc + 1, self.nspecies**3, 5)

        for (in3, num3), (in2, num2), (in1, num1) in product(
                enumerate(infodat.numbers),
                enumerate(infodat.numbers),
                enumerate(infodat.numbers)):
            idx = self.three_center.icon[in1, in2, in3] - 1
            nnz = self.three_center.im[in1, in2]
            for isorp, it in product(
                    range(1, infodat.numshells[num3] + 1),
                    range(5)):
                ffname = f"den3_{it:02}_{isorp:02}"
                ffname += f".{num1:02}.{num2:02}.{num3:02}.dat"
                data = load_data(join(fpath, ffname), header=11)
                ym, ny = read_line(data, "float", "int")
                xm, nx = read_line(data, "float", "int")
                skip_lines(data, lines=5)

                mden[1, isorp, idx], nden[1, isorp, idx] = ym, ny
                mden[0, isorp, idx], nden[0, isorp, idx] = xm, nx
                for iy in range(ny):
                    den[:nnz, :nx, iy, it, isorp, idx] = \
                        read_float_array(data, lines=nx).transpose()[:nnz, :]

        self.three_center.nden = nden
        self.three_center.mden = mden
        self.three_center.den = den

        # Interaction 4
        dens = alloc_float(self.nxm, self.nym, self.me3c_max,
                           self.isorpmax_xc + 1, self.nspecies**3, 5)

        for (in3, num3), (in2, num2), (in1, num1) in product(
                enumerate(infodat.numbers),
                enumerate(infodat.numbers),
                enumerate(infodat.numbers)):
            idx = self.three_center.icon[in1, in2, in3] - 1
            nnz = self.three_center.ims[in1, in2]
            for isorp, it in product(
                    range(1, infodat.numshells[num3] + 1),
                    range(5)):
                ffname = f"deS3_{it:02}_{isorp:02}"
                ffname += f".{num1:02}.{num2:02}.{num3:02}.dat"
                data = load_data(join(fpath, ffname), header=11)
                ym, ny = read_line(data, "float", "int")
                xm, nx = read_line(data, "float", "int")
                skip_lines(data, lines=5)

                for iy in range(ny):
                    dens[:nnz, :nx, iy, it, isorp, idx] = \
                        read_float_array(data, lines=nx).transpose()[:nnz, :]

        self.three_center.dens = dens
