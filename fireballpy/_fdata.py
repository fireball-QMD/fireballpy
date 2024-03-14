from dataclasses import dataclass, field
from os.path import join
import numpy as np

from infodat import InfoDat


@dataclass
class _one_center:
    exc1c0: np.ndarray[int, float] = field(init=False)
    nuxc1c: np.ndarray[int, float] = field(init=False)
    dexc1c: np.ndarray[int, float] = field(init=False)
    dnuxc1c: np.ndarray[int, float] = field(init=False)
    # d2exc1c: np.ndarray[int, float] = field(init=False)
    # d2nuxc1c: np.ndarray[int, float] = field(init=False)


@dataclass
class _fdata:
    nsh_max: int
    nshPP_max: int
    isorpmax: int
    isorpmax_xc: int
    nspecies: int

    one_center: _one_center = field(init=False)

    def load_1c(self, infodat: InfoDat, fpath: str) -> None:
        # Init arrays
        exc1c0 = np.zeros(
            (self.nspecies, self.nsh_max, self.nsh_max),
            dtype=float
        )
        nuxc1c = np.zeros(
            (self.nspecies, self.nsh_max, self.nsh_max),
            dtype=float
        )
        dexc1c = np.zeros(
            (self.nspecies, self.nsh_max, self.nsh_max, self.nsh_max),
            dtype=float
        )
        dnuxc1c = np.zeros(
            (self.nspecies, self.nsh_max, self.nsh_max, self.nsh_max),
            dtype=float
        )
        # d2exc1c = np.zeros(
        #     (self.nspecies, self.nsh_max, self.nsh_max),
        #     dtype=float
        # )
        # d2nuxc1c = np.zeros(
        #     (self.nspecies, self.nsh_max, self.nsh_max,
        #     self.nsh_max, self.nsh_max),
        #     dtype=float
        # )

        for i, num in enumerate(infodat.numbers):
            with open(join(fpath, f"xc1c_dq1.{num:02}.dat"), "r") as fp:
                raw1 = fp.read().splitlines()
            with open(join(fpath, f"exc1crho.{num:02}.dat"), "r") as fp:
                raw2 = fp.read().splitlines()
            with open(join(fpath, f"nuxc1crho.{num:02}.dat"), "r") as fp:
                raw3 = fp.read().splitlines()

            # Get rid of messages
            raw1 = raw1[6:]
            raw2 = raw2[6:]
            raw3 = raw3[6:]

            numsh1 = int(raw1[0].split()[1])
            for j in range(numsh1):
                exc1c0[i, j, :] = list(
                    map(float, raw1[j+1].replace("D", "E").split()))
                nuxc1c[i, j, :] = list(
                    map(float, raw1[j+numsh1+2].replace("D", "E").split()))

            numsh2 = 0
            numsh3 = 0
            for j in range(infodat.numshells[num]):
                tsh2 = j*(numsh2 + 1)
                tsh3 = j*(numsh3 + 1)
                numsh2 = int(raw2[0 + tsh2].split()[1])
                numsh3 = int(raw3[0 + tsh3].split()[1])
                for k in range(numsh2):
                    dexc1c[i, k, :, j] = map(float, raw2[1+k+tsh2])
                for k in range(numsh3):
                    dnuxc1c[i, k, :, j] = map(float, raw3[1+k+tsh3])

        self.one_center = _one_center(
            exc1c0=exc1c0,
            nuxc1c=nuxc1c,
            dexc1c=dexc1c,
            dnuxc1c=dnuxc1c,
        )
