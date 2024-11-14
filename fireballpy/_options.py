from __future__ import annotations
import numpy as np

_icharge_table = {'lowdin': np.int64(1), 'mulliken': np.int64(2),
                 'weighted_lowdin': np.int64(3), 'mulliken_dipole': np.int64(4),
                 'mulliken_dipole_preserving': np.int64(7)}

_idipole_table = {'improved': np.int64(1), 'legacy': np.int64(0)}

_imixer_table = {'anderson': np.int64(1), 'johnson': np.int64(2)}

def get_icharge(charge_method: str) -> np.int64:
    if charge_method not in _icharge_table:
        raise ValueError("Parameter 'charge_method' must be one of "
                         f"{', '.join(_icharge_table)}. Got {charge_method}")
    return _icharge_table[charge_method]


def get_idipole(dipole_method: str) -> np.int64:
    if dipole_method not in _idipole_table:
        raise ValueError("Parameter 'dipole_method' must be one of "
                         f"{', '.join(_idipole_table)}. Got {dipole_method}")
    return _idipole_table[dipole_method]


def get_imixer(mixer_method: str) -> np.int64:
    if mixer_method not in _imixer_table:
        raise ValueError("Parameter 'mixer[\"method\"]' must be one of "
                         f"{', '.join(_imixer_table)}. Got {mixer_method}")
    return _imixer_table[mixer_method]
