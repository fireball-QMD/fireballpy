from typing import Dict, List, Self
from dataclasses import dataclass, field
from os.path import splitext
import json
import tomllib

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
    cutoffs: Dict[int, List[float]]
    qneutral: Dict[int, List[float]]
    wffiles: Dict[int, List[str]]
    nafiles: Dict[int, List[str]]
    energy: Dict[int, float]

    nspecies: int = field(init=False)
    numbers: List[int] = field(init=False)
    numshells: Dict[int, int] = field(init=False)
    numshellsPP: Dict[int, int] = field(init=False)
    elements: Dict[int, Element] = field(init=False)

    def __post_init__(self) -> None:
        if not (self.shells.keys() ==
                self.shellsPP.keys() ==
                self.cutoffPP.keys() ==
                self.cutoffs.keys() ==
                self.qneutral.keys() ==
                self.wffiles.keys() ==
                self.nafiles.keys() ==
                self.energy.keys()):
            raise ValueError("The keys of the dictionaries are not equal")
        numbers = list(self.shells.keys())
        numbers.sort()
        self.numbers = numbers
        self.nspecies = len(numbers)

        numshells = {}
        numshellsPP = {}
        for num in self.numbers:
            if not (len(self.shells[num]) ==
                    len(self.cutoffs[num]) ==
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
                f"Specie {i+1}/{len(self.numbers)} " +
                f"({self.elements[num].name}):"
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
            string.append(
                "- Radial cutoffs (atomic units): " +
                ", ".join(f"{x}" for x in self.cutoffs[num])
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

    @classmethod
    def load(cls, filename: str) -> Self:
        """Load the information of a InfoDat from a .json file or a .toml file.

        Parameters
        ----------
        filename : str
            Name of the file from where InfoDat will be loaded. It must be
            either a .json file or a .toml file

        Raises
        ------
        ValueError
            If the file extension is neither .json nor .toml
        """
        ext = splitext(filename)[1]
        if ext == ".json":
            with open(filename, "w") as fp:
                infodict = json.load(fp)
        elif ext == ".toml":
            with open(filename, "w") as fp:
                infodict = tomllib.load(fp)
        else:
            raise ValueError(f"File extension expected to be either .json or \
                    .toml. Found {ext}")
        return cls(infodict["shells"],
                   infodict["shellsPP"],
                   infodict["cutoffPP"],
                   infodict["cutoffs"],
                   infodict["qneutral"],
                   infodict["wffiles"],
                   infodict["nafiles"],
                   infodict["energy"])

    def save(self, filename: str) -> None:
        """Save the information of this InfoDat into a .json file or
        a .toml file.

        In general a JSON file is recomended for solely storage purposes for
        its cross-compatibility and fast parse. On the other hand, TOML is
        highly recomended if some manual modification is planned

        Parameters
        ----------
        filename : str
            Name of the file where InfoDat will be saved. It must be either a
            .json file or a .toml file

        Raises
        ------
        ValueError
            If the file extension is neither .json nor .toml
        """

        infodict = {
            "shells": self.shells,
            "shellsPP": self.shellsPP,
            "cutoffPP": self.cutoffPP,
            "cutoffs": self.cutoffs,
            "qneutral": self.qneutral,
            "wffiles": self.wffiles,
            "nafiles": self.nafiles,
            "energy": self.energy
        }
        ext = splitext(filename)[1]
        if ext == ".json":
            with open(filename, "w") as fp:
                json.dump(infodict, fp)
        elif ext == ".toml":
            with open(filename, "w") as fp:
                tomllib.dump(infodict, fp)
        else:
            raise ValueError(f"File extension expected to be either .json or \
                    .toml. Found {ext}")


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
    1:  0.2,
    6:  0.87,
    7:  0.78,
    8:  0.71,
    16: 1.24
}

default_cutoffs = {
    1:  [1.0, 0.0],
    6:  [2.0, 2.0],
    7:  [2.0, 3.0, 0.0],
    8:  [2.0, 4.0, 0.0],
    16: [2.0, 4.0, 0.0]
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
    1:  0.0,
    6:  0.0,
    7:  0.0,
    8:  0.0,
    16: 0.0
}

default_infodat = InfoDat(default_shells,
                          default_shellsPP,
                          default_cutoffsPP,
                          default_cutoffs,
                          default_qneutral,
                          default_wffiles,
                          default_nafiles,
                          default_energy)
