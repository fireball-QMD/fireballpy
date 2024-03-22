from dataclasses import dataclass, field
from itertools import product
from os.path import join
import numpy as np

from infodat import InfoDat
from utils import (
    alloc_float,
    alloc_int,
    load_data,
    read_float_array,
    read_int_entry,
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
    exc1: np.ndarray[int, float]
    nuxc: np.ndarray[int, float]
    dexc: np.ndarray[int, float]
    dnuxc: np.ndarray[int, float]


@dataclass
class _TwoCenter:
    numz: np.ndarray[int, int] = field(init=False)
    cl_pp: np.ndarray[int, float] = field(init=False)
    zmax: np.ndarray[int, float] = field(init=False)
    xint: np.ndarray[int, float] = field(init=False)


@dataclass
class _FData:
    nsh_max: int
    nshPP_max: int
    isorpmax: int
    isorpmax_xc: int
    nspecies: int

    one_center: _OneCenter = field(init=False)
    two_center: _TwoCenter = field(init=False)

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
            _ = data.popleft()
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

        self.one_center = _OneCenter(
            exc=exc,
            nuxc=nuxc,
            dexc=dexc,
            dnuxc=dnuxc,
        )

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
                        zm, nz = read_float_array(data)
                        nz = int(round(nz))
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
                            read_float_array(data, lines=nz).T

            self.two_center.numz = numz
            self.two_center.cl_pp = cl_pp
            self.two_center.zmax = zmax
            self.two_center.xint = xint
