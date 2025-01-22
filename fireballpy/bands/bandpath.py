from __future__ import annotations
from copy import deepcopy
from typing import SupportsInt
import re
import numpy as np
from numpy.linalg import norm
from numpy.typing import ArrayLike
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from fireballpy._errors import type_check
from fireballpy.atoms import AtomSystem
from fireballpy.kpoints import KPoints
from fireballpy.fireball import BaseFireball, fbobj_from_obj

from ase import Atoms


class BandPath:
    """Object to create a bandpath plot from a Fireball DFT calculation.

    The implementation is based on ASE's :class:`ase.dft.kpoints.BandPath`
    but it is specialized for Fireball computations, allowing for faster
    and simpler bandpath computation.
    As such, it will only accept objects which have been interfaced with
    fireballpy.
    The usage of this class will be the same for all objects interfaced
    with Fireball.

    Parameters
    ----------
    path : str | ArrayLike[float]
        This parameter may be a ``str`` with the points in the path in which case
        we will try to assign the corresponding k-points.
        The string should look as the collection of all the point names in the path
        in order, using commas to separate non-contiguous paths.
        For more details on how each point may be identified see the documentation
        for the ``special_points`` parameter.
        It may also be a nkpts x 3 array with the points to consider in the path.
        Note that if specified in this form, the path must be continuous.
    npoints : int
        Number of k-points that will be computed for the path. This number must be
        equal or bigger than the number of k-points specified in the ``path`` property.
    atoms : ase.Atoms | None
        An :class:`ase.Atoms` object with a :class:`Fireball` attached as calculator.
        The SCF loop must be computed before.
        This parameter is not compatible with ``fbobj``.
    fbobj : BaseFireball | None
        Directly a low-level :class:`BaseFireball` object.
        In general this is reserved only for advanced users.
        The SCF loop must be computed before.
        This parameter is not compatible with ``atoms``.
    special_points : dict | None:
        Dictionary with the map from the number of a point to its coordinates in
        units of the reciprocal cell.
        Each of the points is defined by an uppercase letter, an arbitrary number of
        lowercase letters and an arbitrary number of numbers.
        If this parameter is not provided, then we will rely on
        :func:`ase.dft.kpoints.get_special_points` to identify the form of the cell
        and the mapping of point names to their corresponding coordinates.

    Methods
    -------
    plot()
        Plot the bandpath using :mod:`matplotlib`.

    Raises
    ------
    ValueError
        If more than one fireballpy object is passed.
    RuntimeError
        If the SCF loop was not computed.
    """
    def __init__(self, *,
                 path: str | ArrayLike,
                 npoints: SupportsInt,
                 atoms: Atoms | None = None,
                 fbobj: BaseFireball | None = None,
                 special_points: dict | None = None) -> None:

        fbobj = fbobj_from_obj(fbobj=fbobj, atoms=atoms)
        if not hasattr(fbobj, 'run_scf') or not fbobj.scf_computed:
            raise RuntimeError("SCF loop needs to be carried first")

        type_check(npoints, SupportsInt, 'npoints')
        self.npoints = int(npoints)

        # See if we have a path as str or defined as an array
        if isinstance(path, str):
            if special_points is None:
                # Get them from ASE
                from ase.dft.kpoints import get_special_points
                self.special_points = get_special_points(fbobj.atomsystem.cell)
            else:
                type_check(special_points, dict, 'special_points')
                self.special_points = special_points

            # Now parse the path
            self.paths = [[name for name in re.split(r'([A-Z][a-z0-9]*)', p) if name] for p in path.split(',')]
        elif isinstance(path, (list, tuple, np.ndarray)):
            # Coordinates
            path = np.ascontiguousarray(path, dtype=np.float64)
            if path.shape != (path.shape[0], 3):
                raise ValueError("Provided ``path`` as coordinates are not a nkpts x 3 array.")

            self.special_points = {f"Kpt{i}": k for i, k in enumerate(path)}
            self.paths = [[k for k in self.special_points]]
        else:
            raise TypeError("Parameter ``path`` could not be identified.")

        for p in self.paths:
            if len(p) < 2:
                raise ValueError("Paths must have at least 2 elements.")

        self.pathcoords = [[self.special_points[p] for p in pp] for pp in self.paths]
        self._paths2kpts(fbobj.atomsystem)

        # Compute new eigenvalues
        self.fermi = fbobj.fermi_level
        self.fbobj = BaseFireball(atomsystem=fbobj.atomsystem, fdatafiles=fbobj.fdatafiles,
                                  kpoints=KPoints(kpts=self.pathcoords, atomsystem=fbobj.atomsystem),
                                  total_charge=fbobj.total_charge, correction=fbobj.correction,
                                  charges_method=fbobj.charges_method, dipole_method=fbobj.dipole_method,
                                  initial_charges=fbobj.shell_charges, mixer_kws=fbobj.mixer_kws)
        self.fbobj.run_scf(fix_charges=True)
        self.eigenvalues = self.fbobj.eigenvalues
        self.nkpts, self.nbands = self.eigenvalues.shape

    # Reference: https://gitlab.com/ase/ase/-/blob/master/ase/dft/kpoints.py?ref_type=heads
    def _paths2kpts(self, atomsystem: AtomSystem) -> None:
        points = np.concatenate(self.pathcoords)
        dists = points[1:] - points[:-1]
        lengths = [norm(d) for d in np.dot(dists, atomsystem.icell)]

        i = 0
        for path in self.pathcoords[:-1]:
            i += len(path)
            lengths[i-1] = np.float64(0.0)
        totallength = sum(lengths)

        i = 0
        x0 = 0.0
        self.pathcoords = []
        self.pathidx = [0]
        for point, dist, length in zip(points[:-1], dists, lengths):
            diff = totallength - x0
            if abs(diff) < 1e-6:
                n = 0
            else:
                n = max(2, int(round(length * (self.npoints - len(self.pathcoords)) / diff)))
            for t in np.linspace(0, 1, n)[:-1]:
                self.pathcoords.append(point + t * dist)
            x0 += length
            if length > 1e-6:
                self.pathidx.append(self.pathidx[-1] + n - 1 + i)
            else:
                i += 1
        if self.pathidx[-1] == self.npoints:
            self.pathidx[-1] -= 1
        self.pathcoords.append(points[-1])
        self.pathcoords = np.ascontiguousarray(self.pathcoords, dtype=np.float64)

    def plot(self, *,
             emin: float | None = None,
             emax: float | None = None,
             bands_fmt: str = '-k',
             fermi_level_fmt: str = '--r',
             center: bool = False,
             ax: Axes | None = None,
             figsize: tuple[float, float] | None = None,  # Parámetro figsize opcional
             show_grid: bool = False,
             line_kws: dict | None = None,
             font_kws: dict | None = None) -> Axes:
        """Plot the bandpath using :mod:`matplotlib`.

        Parameters
        ----------
        emin : float | None
            Lower y-axis bound for the plot.
        emax : float | None
            Upper y-axis bound for the plot.
        bands_fmt : str
            Format string for the band lines. Default is ``'-k'``
            (continous black lines).
        fermi_level_fmt : str
            Format string for the Fermi level. Default is ``'--r'``
            (discontinous red line). To avoid the printing of the
            line use the empty string: ``''``.
        center: bool
            Whether to shift the eigenvalues so that the Fermi level
            lands at the zero energy mark (``True``) or not (``False``, default).
        ax: plt.Axes | None
            If an instance of :class:`matplotlib.pyplot.Axes` is provided, it will
            be used to draw the bandpath on it.
        line_kws: dict | None
            Extra arguments to provide to the calls of :meth:`matplotlib.pyplot.Axes.plot`.
        font_kws: dict | None
            Extra arguments to provide to the font-related methods.
        """

        center = bool(center)
        line_kws = line_kws if line_kws else {}
        font_kws = font_kws if font_kws else {}
        if ax is None:
            default_figsize = (6.9, 5)
            _, ax = plt.subplots(figsize=figsize if figsize else default_figsize, layout='constrained')
            assert ax is not None

        # Define xticks
        xticks = deepcopy(self.paths[0])
        if len(self.paths) > 1:
            for path in self.paths[1:]:
                xticks[-1] += f',{path[0]}'
                xticks.extend(path[1:])
        xticks = list(map(lambda x: x.replace('G', r'$\Gamma$'), xticks))


        # Make the plot
        xplot = np.arange(self.nkpts)
        yplot = self.eigenvalues.T.copy()
        if center:
            yplot -= self.fermi

        for y in yplot:
            ax.plot(xplot, y, bands_fmt, **line_kws)
        if center:
            if fermi_level_fmt:
                ax.plot([xplot[0], xplot[-1]], [0.0, 0.0], fermi_level_fmt, **line_kws)
            ax.set_ylabel('E - E${}_F$ [eV]', **font_kws)
        else:
            if fermi_level_fmt:
                ax.plot([xplot[0], xplot[-1]], [self.fermi, self.fermi], fermi_level_fmt, **line_kws)
            ax.set_ylabel('Energy [eV]', **font_kws)
        ax.set_xlim(xplot[0], xplot[-1])
        ax.set_ylim(bottom=emin, top=emax)
        ax.set_xlabel('')
        ax.set_xticks(self.pathidx, xticks, **font_kws)
        if show_grid:
            grid_style = {'color': 'gray', 'linestyle': '-', 'linewidth': 0.5}
            ax.grid(True, axis='x', **grid_style)  # Mostrar solo líneas verticales (eje x)
        return ax
