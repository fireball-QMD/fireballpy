from abc import ABC, abstractmethod
import numpy as np
from tqdm import tqdm

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from ase import Atoms
from ase.io import read
from ase.units import kcal, mol
from dftd3.ase import DFTD3

from fireballpy import Fireball


KCALMOL = kcal/mol


def parse_selection(sel: str) -> list[int]:
    sel = sel.replace(' ', '')
    pieces = sel.split(',')

    idxs = []
    for piece in pieces:
        if '-' not in piece:
            idxs.append(int(piece) - 1)
            continue
        start, end = piece.split('-')
        for i in range(int(start) - 1, int(end)):
            idxs.append(i)
    return idxs


class BaseScorer(ABC):
    @abstractmethod
    def fdata_score(self, reference: str, fdata_path: str, charges_method: str, verbose: bool = False) -> float:
        ...

    @abstractmethod
    def d3_score(self, reference: str, correction: dict, fdata_path: str, charges_method: str, verbose: bool = False) -> float:
        ...


class InteractionEnergyScorer(BaseScorer):
    _method = 'interaction_energies'

    def __init__(self, data: dict):
        super().__init__()

        self.name = data['name']
        if data['method'].lower() != self._method:
            raise ValueError("Provided dataset is not of the appropiate score type. "
                             f"Got {data['method'].lower()}, expected {self._method}.")

        self.units = 1.0
        if 'units' in data:
            if data['units'].lower() == 'kcal/mol':
                self.units = KCALMOL

        self.curves = data['curves']
        self.fdata_save = {}
        self.d3_save = {}

    def fdata_score(self, reference: str, fdata_path: str, charges_method: str, verbose: bool = False) -> float:
        score = 0.0
        save = {}
        for c, curve in tqdm(self.curves.items(), desc="Computing self-energies", unit='curves', disable=not verbose):
            curvescore = 0.0
            save[c] = {'block_total': 0.0, 'all': [], 'property': []}
            for i, (file, ref) in enumerate(zip(curve['files'], curve['references'][reference])):
                xyz = read(file)
                assert isinstance(xyz, Atoms)

                # First curve: compute energies of each block
                if i == 0:
                    for b, block in curve['blocks'].items():
                        xyz_b = xyz.copy()[parse_selection(block['selection'])]
                        assert isinstance(xyz_b, Atoms)
                        xyz_b.calc = Fireball(fdata='custom', fdata_path=fdata_path,
                                              total_charge=block.get('charge', 0),
                                              charges_method=charges_method,
                                              correction=False, lazy=False)
                        save[c][b] = float(xyz_b.get_potential_energy())
                        save[c]['block_total'] += save[c][b]

                # RMSD for each curve
                xyz.calc = Fireball(fdata='custom', fdata_path=fdata_path,
                                    total_charge=curve.get('charge', 0),
                                    charges_method=charges_method,
                                    correction=False, lazy=False)
                save[c]['all'].append(float(xyz.get_potential_energy()))
                save[c]['property'].append(save[c]['all'][-1] - save[c]['block_total'])
                curvescore += (save[c]['property'][-1] - ref*self.units)**2
            curvescore = np.sqrt(curvescore/len(curve['files']))
            score += curvescore
        score /= len(self.curves)
        self.fdata_save = save
        return score

    def d3_score(self, reference: str, correction: dict, fdata_path: str, charges_method: str, verbose: bool = False) -> float:
        if not self.fdata_save:
            self.fdata_score(reference=reference, fdata_path=fdata_path, charges_method=charges_method)

        score = 0.0
        save = {}
        for c, curve in self.curves.items():
            curvescore = 0.0
            save[c] = {'block_total': 0.0, 'all': [], 'property': []}
            for i, (file, ref) in enumerate(zip(curve['files'], curve['references'][reference])):
                xyz = read(file)
                assert isinstance(xyz, Atoms)

                # First curve: compute energies of each block
                if i == 0:
                    for b, block in curve['blocks'].items():
                        xyz_b = xyz.copy()[parse_selection(block['selection'])]
                        assert isinstance(xyz_b, Atoms)
                        xyz_b.calc = DFTD3(**correction)
                        save[c][b] = float(xyz_b.get_potential_energy()) + self.fdata_save[c][b]
                        save[c]['block_total'] += save[c][b]

                # RMSD for each curve
                xyz.calc = DFTD3(**correction)
                save[c]['all'].append(float(xyz.get_potential_energy()) + self.fdata_save[c]['all'][i])
                save[c]['property'].append(save[c]['all'][-1] - save[c]['block_total'])
                curvescore += (save[c]['property'][-1] - ref*self.units)**2
            curvescore = np.sqrt(curvescore/len(curve['files']))
            score += curvescore
        score /= len(self.curves)
        self.d3_save = save
        return score

    def plot(self, reference: str, nrows: int = 4, ncols: int = 4) -> tuple[list[Figure], list[Axes]]:
        if not self.fdata_save:
            raise RuntimeError("Score not computed previously.")

        ncurves = len(self.curves)
        axsxfig = nrows*ncols
        nfigs = ncurves//axsxfig + 1*(ncurves%axsxfig != 0)

        figs, axss = [], []
        for i, (c, curve) in enumerate(self.curves.items()):
            subpos = i%axsxfig
            if subpos == 0:
                figs.append(plt.figure(figsize=(10, 8), layout='constrained'))
                axss.append([])

                title = self.name
                if nfigs > 1:
                    title += f' ({i+1}:{min((i+1)*axsxfig, ncurves)})'
                title += ': ○ ' + reference + '; □ FB'
                if self.d3_save:
                    title += '; △ FB+D3'
                figs[-1].suptitle(title)

            ax = figs[-1].add_subplot(nrows, ncols, subpos + 1)
            ax.plot(curve['x'], list(map(lambda x: x*self.units, curve['references'][reference])), color='#1E88E5', marker='o', markersize=10)
            ax.plot(curve['x'], self.fdata_save[c]['property'], color='#D81B60', marker='s', markersize=10)
            if self.d3_save:
                ax.plot(curve['x'], self.d3_save[c]['property'], color='#004D40', marker='^', markersize=10)
            ax.set_title(c)
            axss[-1].append(ax)

        return figs, axss
