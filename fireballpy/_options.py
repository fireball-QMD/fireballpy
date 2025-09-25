from __future__ import annotations

_icharge_table = {'lowdin': 1, 'mulliken': 2,
                 'weighted_lowdin': 3, 'mulliken_dipole': 4, 'stationary_charges': 6,
                 'mulliken_dipole_preserving': 7}

_idipole_table = {'improved': 1, 'legacy': 0}

_imixer_table = {'anderson': 1, 'johnson': 2}

def get_icharge(charge_method: str) -> int:
    if charge_method not in _icharge_table:
        raise ValueError("Parameter 'charge_method' must be one of "
                         f"{', '.join(_icharge_table)}. Got {charge_method}")
    return _icharge_table[charge_method]


def get_idipole(dipole_method: str) -> int:
    if dipole_method not in _idipole_table:
        raise ValueError("Parameter 'dipole_method' must be one of "
                         f"{', '.join(_idipole_table)}. Got {dipole_method}")
    return _idipole_table[dipole_method]


def get_imixer(mixer_method: str) -> int:
    if mixer_method not in _imixer_table:
        raise ValueError("Parameter 'mixer[\"method\"]' must be one of "
                         f"{', '.join(_imixer_table)}. Got {mixer_method}")
    return _imixer_table[mixer_method]
