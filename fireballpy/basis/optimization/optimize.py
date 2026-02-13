from abc import ABC, abstractmethod
import json

import optuna
import matplotlib.pyplot as plt

from score import InteractionEnergyScorer


scorers = {
    'interaction_energies': InteractionEnergyScorer,
}


class Optimizer(ABC):
    def __init__(self, datasets: str | list[str], references: str | list[str],
                 charges_method: str, weights: list[float] | None = None) -> None:
        if isinstance(datasets, str):
            if not isinstance(references, str):
                raise ValueError("If ``datasets`` is a ``str`` then ``references`` needs to be a ```str``")
            datasets = [datasets]
            references = [references]
            weights = [1.0]
        elif isinstance(datasets, list):
            if not weights:
                weights = [1.0 for _ in range(len(datasets))]
            if len(weights) != len(datasets):
                raise ValueError("There must be the same elements of ``datasets`` as ``weights``.")
            if len(references) != len(datasets):
                raise ValueError("There must be the same elements of ``datasets`` as ``references``.")
        else:
            raise ValueError("``datasets`` must be either a ``str`` or a list of ``str``")

        self.datasets = datasets
        self.references = references
        self.weights = weights
        self.charges_method = charges_method

        self.scorers = []
        for dataset in self.datasets:
            with open(dataset, 'r', encoding='utf-8') as fp:
                data = json.load(fp)
            try:
                self.scorers.append(scorers[data['method']](data))
            except KeyError:
                raise ValueError(f"Method '{data['method']}' in {dataset} not recognized.")

    def __getitem__(self, key):
        return self.scorers[key]

    @abstractmethod
    def __call__(self, trial: optuna.Trial) -> float:
        ...


class D3Optimizer(Optimizer):
    def __init__(self, datasets: str | list[str], references: str | list[str],
                 charges_method: str, fdata_path: str,
                 damping: str, params: dict,
                 weights: list[float] | None = None, verbose: bool = False, **kwargs) -> None:
        super().__init__(datasets, references, charges_method, weights)
        self.fdata_path = fdata_path
        self.damping = damping
        self.params = params
        self.verbose = verbose
        self.best_params: dict[str, float] = {}

        # Compute the Fireball once
        for scorer, reference in zip(self.scorers, self.references):
            scorer.fdata_score(reference=reference, fdata_path=self.fdata_path,
                               charges_method=self.charges_method, verbose=self.verbose, **kwargs)

    def __call__(self, trial: optuna.Trial) -> float:
        correction = {'damping': self.damping,
                      'params_tweaks': {k: trial.suggest_float(k, **p) for k, p in self.params.items()}}
        score = 0.0
        for weight, scorer, reference in zip(self.weights, self.scorers, self.references):
            score += weight*scorer.d3_score(reference=reference, correction=correction,
                                            fdata_path=self.fdata_path, charges_method=self.charges_method)
        return score


def optimize_dftd3(*, datasets: str | list[str], references: str | list[str],
                   charges_method: str, fdata_path: str,
                   damping: str, params_tweaks: dict,
                   weights: list[float] | None = None, opt_kwargs: dict | None = None,
                   fb_kwargs: dict | None = None) -> tuple[D3Optimizer, optuna.study.Study]:
    """Optimize the DFT-D3 correction of a custom FData.

    Parameters
    ----------
    datasets : str | list[str]
        File(s) with the appropriate JSON specification (see) with the dataset
        information for the benchmark.
    references : str | list[str]
        Name(s) of the reference value for computing the score.
    charges_method : str
        Selected charge method for the SCF. Will only be used to perform the initial
        Fireball computation.
    fdata_path : str
        Selected FData for the initial Fireball computation.
    damping : str
        Damping method for DFT-D3. For more information, see
        `Simple DFT-D3 documentation <https://dftd3.readthedocs.io>`_.
    params_tweaks : dict
        Dictionary where the keys are the names of the parameters
        (see `Simple DFT-D3 documentation <https://dftd3.readthedocs.io>`_) and
        the values are dictionaries with the parameters that may be accepted by
        :method:`optuna.trial.Trial.suggest_float`.
    opt_kwargs
        All parameters accepted by :method:`optuna.study.Study.optimize`
    fb_kwargs
        All parameters accepted by Fireball
    """

    if not opt_kwargs:
        opt_kwargs = {}
    if not fb_kwargs:
        fb_kwargs = {}
    opt = D3Optimizer(datasets=datasets, references=references, charges_method=charges_method,
                      fdata_path=fdata_path, damping=damping, params=params_tweaks,
                      weights=weights,
                      verbose=opt_kwargs.get('show_progress_bar', False), **fb_kwargs)
    study = optuna.create_study()
    study.optimize(opt, **opt_kwargs)
    opt.best_params = study.best_params
    return opt, study


if __name__ == '__main__':
    opt, std = optimize_dftd3(datasets='s66x8.json', references='reference_value',
                              charges_method='weighted_lowdin', fdata_path='/home/roldanc/.cache/fireball/fdatas/biology/',
                              damping='d3bj',
                              params_tweaks={'s8': {'low': 0.1, 'high': 5.0, 'step': 0.0001},
                                             'a1': {'low': 0.1, 'high': 5.0, 'step': 0.0001},
                                             'a2': {'low': 0.1, 'high': 5.0, 'step': 0.0001}},
                              opt_kwargs={'n_trials': 5000, 'n_jobs': -1, 'show_progress_bar': True},
                              fb_kwargs={'mixer_kws': {'method': 'johnson', 'mix_order': 3, 'beta': 0.02}})
    print("Best parameters are:", opt.best_params)
    opt[0].plot(reference='reference_value')
    plt.show()
