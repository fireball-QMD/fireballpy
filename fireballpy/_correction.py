import dftd3.interface as d3i

from fireballpy._errors import type_check
from fireballpy.atoms import AtomSystem

# Could be replaced by constants in future versions if
# ase dependency is to be removed
from ase.units import Bohr, Hartree


# Reference: https://stackoverflow.com/a/60769071
class Correction:
    """Low-level interaction with DFT-D3 or DFT-D4 correction

    Parameters
    ----------
    kind : str
        Right now it can only be ``'dftd3'``.
    atomsystem : AtomSystem
        An AtomSystem object with the information of the atomic numbers, positions and the unit cell.
    damping : str | None
        Damping model used. Compulsory for DFT-D3.
    method : str | None
        Named set of parameters to be used. Incompatible with ``params_tweaks``.
    params_tweaks : dict | None
        Dictionary with the explicit values for the parameters. Incompatible with ``method``.
    """
    def __new__(cls, *,
                kind: str,
                atomsystem: AtomSystem,
                damping: str | None = None,
                method: str | None = None,
                params_tweaks: dict | None = None):
        type_check(kind, str, 'kind')
        type_check(atomsystem, AtomSystem, 'atomsystem')
        if kind == 'dftd3':
            type_check(damping, str, 'damping', ' for DFT-D3 correction')
        else:
            raise ValueError("Parameter ``kind`` must be ``'dftd3'``.")
        if method is not None and params_tweaks is not None:
            raise ValueError("Parameters ``method`` and ``params_tweaks`` cannot be set at the same time.")
        if method is None and params_tweaks is None:
            raise ValueError("Parameters ``method`` and ``params_tweaks`` cannot be ``None`` at the same time.")
        if method is None:
            type_check(params_tweaks, dict, 'params_tweaks')
        if params_tweaks is None:
            type_check(method, str, 'method')

        subclass_map = {subclass.kind: subclass for subclass in cls.__subclasses__()}  # type: ignore
        subclass = subclass_map[kind]
        return super(cls, subclass).__new__(subclass)

    def update_coords(self, atomsystem: AtomSystem) -> None:
        self.disp.update(atomsystem.positions/Bohr)  # type: ignore

    def correct(self) -> None:
        self.res = self.disp.get_dispersion(param=self.dpar, grad=True)  # type: ignore
        self.res['energy'] *= Hartree
        self.res['gradient'] *= Hartree/Bohr


class _DFTD3Correction(Correction):
    """Low-level interaction with DFT-D3 correction

    Parameters
    ----------
    atomsystem : AtomSystem
        An AtomSystem object with the information of the atomic numbers, positions and the unit cell.
    damping : str
        Damping model used.
    method : str | None
        Named set of parameters to be used. Incompatible with ``params_tweaks``.
    params_tweaks : dict | None
        Dictionary with the explicit values for the parameters. Incompatible with ``method``.
    """

    kind = 'dftd3'

    _damping_param = {
        "d3bj": d3i.RationalDampingParam,
        "d3zero": d3i.ZeroDampingParam,
        "d3bjm": d3i.ModifiedRationalDampingParam,
        "d3mbj": d3i.ModifiedRationalDampingParam,
        "d3zerom": d3i.ModifiedZeroDampingParam,
        "d3mzero": d3i.ModifiedZeroDampingParam,
        "d3op": d3i.OptimizedPowerDampingParam,
    }

    def __init__(self, kind: str, atomsystem: AtomSystem, damping: str,
                 method: str | None = None, params_tweaks: dict | None = None) -> None:
        self.damping = damping
        self.method = method
        self.params_tweaks = params_tweaks
        self._create_api_calculator(atomsystem)
        self._create_damping_param()

    def _create_api_calculator(self, atomsystem: AtomSystem) -> None:
        """Create a new API calculator object.
        """
        try:
            self.disp = d3i.DispersionModel(atomsystem.numbers,
                                            atomsystem.positions/Bohr)
        except RuntimeError:
            raise ValueError("Cannot construct dispersion model for dftd3")

    def _create_damping_param(self) -> None:
            """Create a new API damping parameter object.
            """
            try:
                params_tweaks = self.params_tweaks if self.params_tweaks is not None else {'method': self.method}
                self.dpar = self._damping_param[self.damping](**params_tweaks)
            except RuntimeError:
                raise ValueError("Cannot construct damping parameter for dftd3")
