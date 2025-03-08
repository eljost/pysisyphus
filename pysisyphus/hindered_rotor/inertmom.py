# [1] https://doi.org/10.1063/1.473958
#     Ab initio statistical thermodynamical models for the
#     computation of third-law entropies
#     East, Radom, 1997

from typing import Literal, Sequence

import numpy as np

from pysisyphus.hessian_proj import get_standard_orientation


def get_com(coords3d: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Get center-of-mass.

    Parameters
    ----------
    coords3d
        2d array of shape (natoms, 3) holding Cartesian atomic coordinates.
    masses
        1d array of shape (natoms, ) holding atomic mass.

    Returns
    -------
    com
        1d array of shape (3, ) holding the center of mass in units of coords3d.
    """
    total_mass = masses.sum()
    com = 1 / total_mass * np.sum(coords3d * masses[:, None], axis=0)
    return com


def inertia_1n(
    coords3d: np.ndarray,
    masses: np.ndarray,
    axis: np.ndarray,
    pivot: np.ndarray,
) -> float:
    """Moment of interia for rotation of coords3d w/ pivot around axis.

    Parameters
    ----------
    coords3d
        2d array of shape (natoms_top, 3) holding Cartesian atomic coordinates. In
        this function, coords3d should only contain the coordinates of the atoms
        in the selected top.
    masses
        1d array of shape (natoms_top, ) holding atomic mass. Similar to coords3d,
        only masses for the atoms in the selected top should be included.
    axis
        1d array of shape (3, ) with norm 1, holding the axis of rotation.
    pivot
        1d array of shape (3, ) holding the pivot point. The centroid of coords3d
        will be shifted to the pivot.

    Returns
    -------
    I
        Moment of intertia for rotation of coords3d around axis with pivot.
    """
    coords3d = coords3d - pivot[None, :]
    # Project out components along axis of rotation.
    coords3d_proj = (
        coords3d - np.einsum("ij,j->i", coords3d, axis)[:, None] * axis[None, :]
    )
    dists2 = np.linalg.norm(coords3d_proj, axis=1) ** 2
    inertia_moment = (masses * dists2).sum()
    return inertia_moment


def get_top_moment_of_inertia(
    coords3d: np.ndarray,
    masses: np.ndarray,
    top_inds: list[int],
    bond_inds: Sequence[int],
    m: Literal[1, 2, 3] = 2,
    n: Literal[1, 2, 3] = 3,
) -> tuple[float, float]:
    """
    Parameters
    ----------
    coords3d
        2d array of shape (natoms, 3) holding the Cartesian atomic coordinates of
        the molecule.
    masses
        1d array of shape (natoms, ) holding atomic mass. Similar to coords3d,
        only masses for the atoms in the selected top should be included.
    top_inds
    axis
        1d array of shape (3, ) with norm 1, holding the axis of rotation.
    pivot
        1d array of shape (3, ) holding the pivot point. The centroid of coords3d
        will be shifted to the pivot.

    Returns


    Returns
    -------
    top_inertia
    rest_inertia
    """
    std_orient = get_standard_orientation(coords3d, masses)
    # Continue w/ coordinates in standard orientation
    coords3d = std_orient.coords3d

    # Inertia tensor, its eigenvalues and -vectors.
    I = std_orient.inertia_tensor
    Iw, Iv = np.linalg.eigh(I)

    top_c3d = coords3d[top_inds]
    top_masses = masses[top_inds]
    rest_inds = [i for i, _ in enumerate(coords3d) if i not in top_inds]
    rest_c3d = coords3d[rest_inds]
    rest_masses = masses[rest_inds]

    top_pivot, rest_pivot = bond_inds
    if top_pivot not in top_inds:
        rest_pivot, top_pivot = bond_inds
    assert (top_pivot in top_inds) and (rest_pivot in rest_inds)
    bond_vec = coords3d[top_pivot] - coords3d[rest_pivot]
    bond_vec /= np.linalg.norm(bond_vec)

    # Rotation around bond axis
    if n == 1:
        top_pivot = coords3d[top_pivot]
        rest_pivot = coords3d[rest_pivot]
        axis = bond_vec
    # Rotation around bond axis, but centered at center of mass
    elif n == 2:
        top_pivot = get_com(top_c3d, top_masses)
        rest_pivot = get_com(rest_c3d, rest_masses)
        axis = bond_vec
    # Rotation around vector between centers of mass, centered at center of mass
    elif n == 3:
        top_pivot = get_com(top_c3d, top_masses)
        rest_pivot = get_com(rest_c3d, rest_masses)
        axis = top_pivot - rest_pivot

    axis /= np.linalg.norm(axis)
    top_inertia = inertia_1n(top_c3d, top_masses, axis, top_pivot)
    rest_inertia = inertia_1n(rest_c3d, rest_masses, axis, rest_pivot)

    if m == 1:
        # No correction for m=1
        pass
    elif m == 2:
        # Eq. (4) in [1]
        inertia = (top_inertia * rest_inertia) / (top_inertia + rest_inertia)
        top_inertia = inertia
        rest_inertia = inertia
    elif m == 3:
        Iv0, Iv1, Iv2 = Iv.T
        dir_cos0 = axis.dot(Iv0)
        dir_cos1 = axis.dot(Iv1)
        dir_cos2 = axis.dot(Iv2)
        # Eq. (6) in [1]
        factor = dir_cos0**2 / Iw[0] + dir_cos1**2 / Iw[1] + dir_cos2**2 / Iw[2]
        # Eq. (5) in [1]
        top_inertia -= top_inertia**2 * factor
        rest_inertia -= rest_inertia**2 * factor

    return top_inertia, rest_inertia
