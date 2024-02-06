from typing import Dict, List, Self
from dataclasses import dataclass, field
from os.path import splitext
import json

from mendeleev import element
from mendeleev.models import Element

@dataclass
class InfoDat:
    """Class to keep all the info present in the classsical info.dat file,
    needed for the adjusting the different parameters of each specie in the
    computation
    """
    shells: Dict[int, List[int]]
    shellsPP: Dict[int, List[int]]
    cutoffPP: Dict[int, float]
    qneutral: Dict[int, List[float]]
    wffiles: Dict[int, List[str]]
    nafiles: Dict[int, List[str]]
    energy: Dict[int, float]

    numbers: List[int] = field(init=False)
    numshells: Dict[int, int] = field(init=False)
    numshellsPP: Dict[int, int] = field(init=False)
    elements: Dict[int, Element] = field(init=False)

    def __post_init__(self) -> None:
        if not (self.shells.keys() ==
                self.shellsPP.keys() ==
                self.cutoffPP.keys() ==
                self.qneutral.keys() ==
                self.wffiles.keys() ==
                self.nafiles.keys() ==
                self.energy.keys()):
            raise ValueError("The keys of the dictionaries are not equal")
        numbers = list(self.shells.keys())
        numbers.sort()
        self.numbers = numbers

        numshells = {}
        numshellsPP = {}
        for num in self.numbers:
            if not (len(self.shells[num]) ==
                    len(self.qneutral[num]) ==
                    len(self.wffiles[num]) ==
                    len(self.nafiles[num]) - 1):
                raise ValueError(
                        f"The number of items in the dictionaries do not match\
                                for element {num}"
                        )
            numshells[num] = len(self.shells[num])
            numshellsPP[num] = len(self.shellsPP[num])
        self.numshells = numshells
        self.numshellsPP = numshellsPP

        elements = {}
        for num in self.numbers:
            elements[num] = element(num)
        self.elements = elements

    def __str__(self) -> str:
        string = []
        for i, num in enumerate(self.numbers):
            string.append(
                f"Specie {i+1}/{len(self.numbers)} ({self.elements[num].name}):"
            )
            string.append(
                "- Angular momentum of the shells: " +
                ", ".join(f"{x}" for x in self.shells[num])
            )
            string.append(
                "- Angular momentum of the shells (Pseudopotential): " +
                ", ".join(f"{x}" for x in self.shellsPP[num])
            )
            string.append(
                f"- Radial cutoff of the pseudopotential: {self.cutoffPP[num]}"
            )
            string.append(f"- Atomic energy: {self.energy[num]}")
            string.append(
                    "- Wavefunction files: " +
                    ", ".join(self.wffiles[num])
            )
            string.append(
                    "- Potential files: " +
                    ", ".join(self.nafiles[num])
            )
            string.append("")
        return "\n".join(string[:-1])

    def save(self, filename: str) -> None:
        jsondict = {
            "shells": self.shells,
            "shellsPP": self.shellsPP,
            "cutoffPP": self.cutoffPP,
            "qneutral": self.qneutral,
            "wffiles": self.wffiles,
            "nafiles": self.nafiles,
            "energy": self.energy
        }
        if splitext(filename)[1] != ".json":
            filename += ".json"
        with open(filename, "w") as fp:
            json.dump(jsondict, fp)

    @classmethod
    def load(cls, filename: str) -> Self:
        with open(filename) as fp:
            jsondict = json.load(fp)
        return cls(jsondict["shells"],
                   jsondict["shellsPP"],
                   jsondict["cutoffPP"],
                   jsondict["qneutral"],
                   jsondict["wffiles"],
                   jsondict["nafiles"],
                   jsondict["energy"])


default_shells = {
    1:  [0, 0],
    6:  [0, 1],
    7:  [0, 1, 2],
    8:  [0, 1, 2],
    16: [0, 1, 2]
}

default_shellsPP = {
    1:  [0],
    6:  [0, 1],
    7:  [0, 1],
    8:  [0, 1],
    16: [0, 1]
}

default_cutoffsPP = {
    1: 0.2,
    6: 0.87,
    7: 0.78,
    8: 0.71,
    16: 1.24
}

default_qneutral = {
    1:  [4.65, 4.65],
    6:  [5.95, 5.95],
    7:  [5.42, 5.42, 4.62],
    8:  [5.32, 5.32, 5.32],
    16: [7.0, 7.0, 5.5]
}

default_wffiles = {
    1:  ["cinput/001_465.wf1", "cinput/001_465.ewf1"],
    6:  ["cinput/006_595.wf1", "cinput/006_595.wf2"],
    7:  ["cinput/007_542.wf1", "cinput/007_542.wf2", "cinput/007_462.wf3"],
    8:  ["cinput/008_532.wf1", "cinput/008_532.wf2", "cinput/008_532.wf3"],
    16: ["cinput/016_700.wf1", "cinput/016_700.wf2", "cinput/016_550.wf3"]
}

default_nafiles = {
    1:  [
            "cinput/001_465.na0",
            "cinput/001_465.na1",
            "cinput/001_465.ena1"],
    6:  [
            "cinput/006_595.na0",
            "cinput/006_595.na1",
            "cinput/006_595.na2"
        ],
    7:  [
            "cinput/007_542.na0",
            "cinput/007_542.na1",
            "cinput/007_542.na2",
            "cinput/007_462.na3"],
    8:  [
            "cinput/008_532.na0",
            "cinput/008_532.na1",
            "cinput/008_532.na2",
            "cinput/008_532.na3"
        ],
    16: [
            "cinput/016_700.na0",
            "cinput/016_700.na1",
            "cinput/016_700.na2",
            "cinput/016_550.na3"
        ]
}

default_energy = {
    1: 0.0,
    6: 0.0,
    7: 0.0,
    8: 0.0,
    16: 0.0
}

default_infodat = InfoDat(default_shells,
                          default_shellsPP,
                          default_cutoffsPP,
                          default_qneutral,
                          default_wffiles,
                          default_nafiles,
                          default_energy)
