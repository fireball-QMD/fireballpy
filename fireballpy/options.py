from __future__ import annotations
from enum import IntEnum

class ChargesMethod(IntEnum):
    LOWDIN = 1
    MULLIKEN = 2
    WEIGHTED_LOWDIN = 3
    MULLIKEN_DIPOLE = 4
    STATIONARY_CHARGES = 6
    MULLIKEN_DIPOLE_PRESERVING = 7

class DipoleMethod(IntEnum):
    LEGACY = 0
    IMPROVED = 1

class MixerMethod(IntEnum):
    ANDERSON = 1
    JOHNSON = 2

class FixShellsMethod(IntEnum):
    NONE = 0
    ONLYD = 1
    AUTO = 2
