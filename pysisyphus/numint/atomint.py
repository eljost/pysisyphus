# [1] https://doi.org/10.1063/1.5049435
#     An improved molecular partitioning scheme for numerical
#     quadratures in density functional theory
#     Laqua, Kussmann, Ochsenfeld, 2018
# [2] https://doi.org/10.1021/ct050190+
#     Distributed Multipole Analysis: Stability for Large Basis Sets
#     Stone, 2005
# [3] https://doi.org/10.1063/1.1725697
#     Atomic Radii in Crystals
#     Slater, 1964


import collections
from typing import Tuple

# import math

import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.numint import angint, radint


# Atomic radii in Bohr as found Table II in [1]
_RADII = {
    "h": 0.8,
    "he": 0.9,
    "li": 1.8,
    "be": 1.4,
    "b": 1.3,
    "c": 1.1,
    "n": 0.9,
    "o": 0.9,
    "f": 0.9,
    "ne": 0.9,
    "na": 1.4,
    "mg": 1.3,
    "al": 1.3,
    "si": 1.2,
    "p": 1.1,
    "s": 1.0,
    "cl": 1.0,
    "ar": 1.0,
    "k": 1.5,
    "ca": 1.4,
    "ga": 1.1,
    "ge": 1.0,
    "as": 0.9,
    "se": 0.9,
    "br": 0.9,
    "kr": 0.9,
    "sc": 1.3,
    "ti": 1.2,
    "v": 1.2,
    "cr": 1.2,
    "mn": 1.2,
    "fe": 1.2,
    "co": 1.2,
    "ni": 1.1,
    "cu": 1.1,
    "zn": 1.1,
}
# Return 1.0 Bohr for all elements heavier than Zn
ATOMIC_RADII = collections.defaultdict(lambda: 1.0, **_RADII)

# Slater-Bragg atomic radii in Bohr as utilized in Stone's GDMA code.
# H uses 0.65 Angstrom. The number below are in Angstrom though but are
# converted to Bohr. Disregarding the value for Hydrogen which was
# originally is 0.25 Angstrom all values are from [3].
_GDMA_RADII = {
    elem: radius * ANG2BOHR
    for elem, radius in {
        "h": 0.65,
        "he": 0.5,
        "li": 1.45,
        "be": 1.05,
        "b": 0.85,
        "c": 0.7,
        "n": 0.65,
        "o": 0.6,
        "f": 0.5,
        "ne": 0.5,
        "na": 1.8,
        "mg": 1.5,
        "al": 1.25,
        "si": 1.1,
        "p": 1.0,
        "s": 1.0,
        "cl": 1.0,
        "ar": 1.0,
        "k": 2.2,
        "ca": 1.8,
        "sc": 1.6,
        "ti": 1.4,
        "v": 1.35,
        "cr": 1.4,
        "mn": 1.4,
        "fe": 1.4,
        "co": 1.35,
        "ni": 1.35,
        "cu": 1.35,
        "zn": 1.35,
        "ga": 1.3,
        "ge": 1.25,
        "as": 1.15,
        "se": 1.15,
        "br": 1.15,
        "kr": 1.15,
        "rb": 2.35,
        "sr": 2.0,
        "y": 1.8,
        "zr": 1.55,
        "nb": 1.45,
        "mo": 1.45,
        "tc": 1.35,
        "ru": 1.3,
        "rh": 1.35,
        "pd": 1.4,
        "ag": 1.6,
        "cd": 1.55,
        "in": 1.55,
        "sn": 1.45,
        "sb": 1.45,
        "te": 1.4,
        "i": 1.4,
        "xe": 1.4,
    }.items()
}
# Return 2.5 Bohr for all elements beyond Xe
GDMA_ATOMIC_RADII = collections.defaultdict(lambda: 2.5, **_GDMA_RADII)


def get_grids_for_fermion_kind(kind: str) -> Tuple[int, int, int, int]:
    """Get number of radial and angular grid points for given kind.

    As outlined in Table III. in [1].

    Parameters
    ----------
    kind
        One of 'g1' to 'g7'.

    Returns
    -------
    npoints
        Tuple of 4 integers comprising the number of radial grid points
        and three angular (Lebedev) grid sizes for the inner, middle and
        outer parts of an atom.
    """
    return {
        # grid: (n_rad, inner, medium, outer)
        # n_rad: (base) number of radial grid points
        # inner/medium/outer: angular grids points
        "g1": (35, 14, 50, 110),
        "g2": (40, 26, 74, 194),
        "g3": (50, 38, 110, 302),
        "g4": (55, 50, 194, 434),
        "g5": (60, 50, 194, 590),
        "g6": (70, 86, 302, 974),
        "g7": (80, 110, 434, 1454),
    }[kind.lower()]


def get_extra_n_rad(atom: str) -> int:
    """Get additional number of radial points for the given atom.

    As outlined (A6) in [1].

    Parameters
    ----------
    atom
        Atomic symbol.

    Returns
    -------
    nextra
        Number of extra radial points for the given atom.
    """
    atomic_number = ATOMIC_NUMBERS[atom.lower()]
    # No extra points for hydrogen & Helium
    if 1 <= atomic_number <= 2:
        n_rad_extra = 0
    # Lithium to Neon
    elif 3 <= atomic_number <= 10:
        n_rad_extra = 5
    # Sodium to Argon
    elif 11 <= atomic_number <= 18:
        n_rad_extra = 10
    # Potassium to Krypton
    elif 19 <= atomic_number <= 36:
        n_rad_extra = 20
    # Rubidium to Xenon
    elif 37 <= atomic_number <= 54:
        n_rad_extra = 25
    # From Caesium onwards
    else:
        n_rad_extra = 30
    return n_rad_extra


def combine_grids(
    origin: np.ndarray, rad_grid: np.ndarray, ang_grid: np.ndarray
) -> np.ndarray:
    """Form atomic grid from cartesian product of radial and angular grid.

    Parameters
    ----------
    origin
        Origin of the grid. Usually the coordinates of the parent atom.
    rad_grid
        Radial grid. 2d-array of shape (nrad, nrad) with the first column
        containing radii and the second column the associated weights.
    ang_grid
        Cartesian angular grid of shape (nang, 4). The first three column
        contain the Cartesian coordiantes (x, y, z) of the grid points; the
        last column contains the associated integration weights.

    Returns
    -------
    grid
        Cartesian grid of shape (nrad * nang, 4). The first three column
        contain the Cartesian coordiantes (x, y, z) of the grid points; the
        last column contains the associated integration weights.
    """
    x0, y0, z0 = origin
    rr, rw = rad_grid.T
    ax, ay, az, aw = ang_grid.T
    grid = np.empty((rw.size * aw.size, 4))
    # XYZ columns
    grid[:, 0] = (x0 + ax[:, None] * rr).flatten()
    grid[:, 1] = (y0 + ay[:, None] * rr).flatten()
    grid[:, 2] = (z0 + az[:, None] * rr).flatten()
    # Weight column
    grid[:, 3] = (aw[:, None] * rw).flatten()
    return grid


def get_fermion_atomic_grid(
    atom: str, origin: np.ndarray, kind: str = "g3"
) -> Tuple[np.ndarray, np.ndarray]:
    """Get atomic grid, as described by the Fermion++ developers.

    Parameters
    ----------
    atom
        Chemical symbol of the atom.
    origin
        Coordinates of the grid's center/origin. Usually the coordiantes
        of the host atom.
    kind
        One of 'g1' to 'g7'. Higher kinds correspond to finer grids.

    Returns
    -------
    xyz
        2d array containing Cartesian gridpoints of shape (npoints, 3).
    weights
        1d array containing integration weights of shape (npoints, ).
    """
    atom = atom.lower()
    atomic_radius = ATOMIC_RADII[atom.lower()]

    # Number of grid points (radial and for three angular grids)
    n_rad, inner, medium, outer = get_grids_for_fermion_kind(kind)

    # Radial grid
    n_rad_extra = get_extra_n_rad(atom)
    n_rad = n_rad + n_rad_extra
    rad_grid = radint.chebyshev_2nd_kind(n_rad, atomic_radius=atomic_radius)
    # Integration in spherical coordinates implies multiplying the integrand by r²
    rad_grid[:, 1] *= rad_grid[:, 0] ** 2

    mask = np.full(n_rad, False, dtype=bool)
    # See eq. (A5) in [1]
    # I certainly don't get the number of grid points per C-atom as given in the
    # last column of Table III in [1] ...
    n_rad_3 = n_rad // 3
    n_rad_2 = n_rad // 2
    # The two lines below result in a lower number of points (~ 2% - 5% less).
    # n_rad_3 = math.ceil(n_rad / 3) + 1
    # n_rad_2 = math.ceil(n_rad / 2) + 1

    inner_mask = mask.copy()
    inner_mask[:n_rad_3] = True
    ninner = inner_mask.sum()
    medium_mask = mask.copy()
    medium_mask[n_rad_3:n_rad_2] = True
    nmedium = medium_mask.sum()
    outer_mask = mask.copy()
    outer_mask[n_rad_2:] = True
    nouter = outer_mask.sum()
    assert sum((ninner, nmedium, nouter)) == n_rad

    inner_ang_grid = angint.get_lebedev_grid(inner)
    inner_grid = combine_grids(origin, rad_grid[inner_mask], inner_ang_grid)
    medium_ang_grid = angint.get_lebedev_grid(medium)
    medium_grid = combine_grids(origin, rad_grid[medium_mask], medium_ang_grid)
    outer_ang_grid = angint.get_lebedev_grid(outer)
    outer_grid = combine_grids(origin, rad_grid[outer_mask], outer_ang_grid)

    grid = np.concatenate((inner_grid, medium_grid, outer_grid), axis=0)
    return grid[:, :3], grid[:, 3]


def get_gdma_atomic_grid(
    origin: np.ndarray, atom: str = "c", n_rad=80, n_ang=590, r_scale=2.0
):
    """Get atomic grid, similar to the original GDMA grid.

    In the original paper [2] Stone proposed 80 radial points w/ Euler-Maclaurin-
    integration and a 590-point Lebedev grid.

    In GDMA the number of expansion sites may/can be different from the atomic sites.
    In such cases we default to the atomic radius of carbon, but the user may choose
    a different atom.

    Parameters
    ----------
    origin
        Coordinates of the grid's center/origin. Usually the coordiantes
        of the host atom.
    atom
        Optional, defaults to carbon.
    n_rad
        Number of radial points, defaults to 80.
    n_ang
        Number of angular points, defaults to 590.
    r_scale
        Scaling factor for the atomic radii.

    Returns
    -------
    xyz
        2d array containing Cartesian gridpoints of shape (npoints, 3).
    weights
        1d array containing integration weights of shape (npoints, ).
    """
    atom = atom.lower()
    atomic_radius = r_scale * GDMA_ATOMIC_RADII[atom.lower()]

    # rad_grid = radint.chebyshev_2nd_kind(n_rad, atomic_radius=atomic_radius)
    rad_grid = radint.euler_maclaurin_dma(n_rad, m_r=2, atomic_radius=atomic_radius)
    # Integration in spherical coordinates implies multiplying the integrand by r²
    rad_grid[:, 1] *= rad_grid[:, 0] ** 2

    ang_grid = angint.get_lebedev_grid(n_ang)

    grid = combine_grids(origin, rad_grid, ang_grid)
    return grid[:, :3], grid[:, 3]


get_atomic_grid = get_fermion_atomic_grid
