from typing import Optional
import dftd3.interface as d3i

from .atoms import AtomSystem

# Could be replaced by constants in future versions if
# ase dependency is to be removed
from ase.units import Bohr

class DFTD3Correction:

    _damping_param = {
        "d3bj": d3i.RationalDampingParam,
        "d3zero": d3i.ZeroDampingParam,
        "d3bjm": d3i.ModifiedRationalDampingParam,
        "d3mbj": d3i.ModifiedRationalDampingParam,
        "d3zerom": d3i.ModifiedZeroDampingParam,
        "d3mzero": d3i.ModifiedZeroDampingParam,
        "d3op": d3i.OptimizedPowerDampingParam,
    }

    def __init__(self,
                 atomsystem: AtomSystem,
                 damping: str,
                 method: Optional[str] = None,
                 params_tweaks: Optional[dict[str, float]] = None) -> None:
        assert isinstance(atomsystem, AtomSystem)
        assert isinstance(damping, str)
        self.atomsystem = atomsystem
        self.damping = damping
        self.method = method
        self.params = params_tweaks
        self._create_api_calculator()
        self._create_damping_param()

    def _create_api_calculator(self) -> None:
        """Create a new API calculator object"""

        try:
            self.disp = d3i.DispersionModel(
                self.atomsystem.nums,
                self.atomsystem.pos/Bohr,
                self.atomsystem.cell/Bohr,
                self.atomsystem.pbc
            )
        except RuntimeError:
            raise ValueError("Cannot construct dispersion model for dftd3")

    def _create_damping_param(self) -> None:
            """Create a new API damping parameter object"""

            if self.method is not None and self.params is not None:
                raise ValueError("'method' and 'params_tweaks' cannot be set at the same time.")
            if self.method is None and self.params is None:
                raise ValueError("'method' and 'params_tweaks' cannot be unset at the same time.")

            try:
                params_tweaks = self.params if self.params is not None else {'method': self.method}
                self.dpar = self._damping_param[self.damping](**params_tweaks)
            except RuntimeError:
                raise ValueError("Cannot construct damping parameter for dftd3")


def new_correction(atomsystem: AtomSystem,
                   corrtype: str,
                   damping: str,
                   method: Optional[str] = None,
                   params_tweaks: Optional[dict[str, float]] = None) -> DFTD3Correction:
    """Wrapper to create the appropiate correction.

    Parameters
    ----------
    atomsystem : AtomSystem
        An AtomSystem object with the information of the atomic numbers, positions and the unit cell.
    corrtype : str
        Right it can only be ``'dftd3'``.
    damping : str
        Damping model used.
    method : Optional[str]
        Named set of parameters to be used. Incompatible with ``params_tweaks``.
    params_tweaks : Optional[dict]
        Dictionary with the explicit values for the parameters. Incompatible with ``method``.

    Returns
    -------
    DFTD3Correction
        Created correction object to be used in Fireball.
    """
    assert isinstance(corrtype, str)
    if method is not None and params_tweaks is not None:
        raise ValueError("'method' and 'params_tweaks' cannot be set at the same time.")
    if method is None and params_tweaks is None:
        raise ValueError("'method' and 'params_tweaks' cannot be unset at the same time.")
    if corrtype not in ['dftd3']:
        raise ValueError("correction type must be 'dftd3'")
    return DFTD3Correction(atomsystem, damping, method, params_tweaks)
