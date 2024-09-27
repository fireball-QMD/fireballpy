from __future__ import annotations
import warnings
from dftd3.ase import DFTD3  # type: ignore

_dftd3_params = {
    'biology': {
        'damping': 'd3bj',
        'params_tweaks': {'s8': 1.80, 'a1': 0.62, 'a2': 2.50}
    },
    'biology_ext': {
        'damping': 'd3bj',
        'params_tweaks': {'s8': 1.90, 'a1': 0.65, 'a2': 3.30}
    }
}


def get_dftd3(name: str) -> DFTD3 | None:
    if name not in _dftd3_params:
        warnings.warn(f'FData {name} does not have optimized parameters '
                      'for DFTD3. No correction will be applied. '
                      'You may consider to add this correction manually. '
                      'Please refer to ``dftd3.ase.DFTD3`` for how to do this.')
        return None
    return DFTD3(**_dftd3_params[name])

