# [1] https://doi.org/10.1016/0009-2614(96)00600-8
#     Achieving linear scaling in exchange-correlation
#     density functional quadratures
#     Stratmann, Scuseria, Frisch, 1996


import numba
import numpy as np

from pysisyphus.numint.atomint import get_atomic_grid


@numba.jit(nopython=True, cache=True)
def stratmann_partitioning(
    coords3d: np.ndarray,
    atomic_grids: numba.typed.List[np.ndarray],
    part_weights: np.ndarray,
    a: float = 0.64,
):
    """Partitioning of atomic weights according to Stratmann, Scuseria and Frisch.

    Parameters
    ----------
    coords3d
        Array holding the coordinates of the grid parents; usually the atomic coordiantes.
    atomic_grid
        Numba typed-List containig atomic grid arrays. The lengths of the respective arrays
        sum to 'npoints'.
    part_weights
        1d array  of shape(npoints, ) containing the partitioning weights that have to be
        multiplied with the integration weights.
    a
        Scaling factor for the smoothing function. Defaults to 0.64.
    """
    assert len(coords3d) == len(atomic_grids)
    ncenters = len(coords3d)
    part_weights[:] = 1.0
    # Return early for a single atom/grid.
    if len(coords3d) == 1:
        return

    # Calculate interatomic distances
    all_center_dists = np.zeros((ncenters, ncenters))
    for i in range(ncenters):
        for j in range(i + 1, ncenters):
            dist = np.linalg.norm(coords3d[i] - coords3d[j])
            all_center_dists[i, j] = dist
            all_center_dists[j, i] = dist

    cur_index = 0
    start_indices = [cur_index]
    for grid in atomic_grids[:-1]:
        cur_index += len(grid)
        start_indices.append(cur_index)

    for p in range(ncenters):
        other_dists = all_center_dists[p]
        # Pick the 2nd lowest distances, as the (first) lowest distance will be the
        # distance of center p to itself (0.0).
        closest_dist = other_dists[other_dists.argsort()[1]]

        # Eq. (15) in [XX]
        thresh = 0.5 * (1.0 - a) * closest_dist
        mus = np.empty((ncenters, ncenters))
        start_index = start_indices[p]
        grid_p = atomic_grids[p]
        for g in range(len(grid_p)):
            # )
            point_g = grid_p[g]
            # Distance of grid point to atomic centers
            dists = np.empty(ncenters)
            for i in range(ncenters):
                dists[i] = np.linalg.norm(coords3d[i] - point_g)
            # Eq. (15) in [1]
            if dists[p] < thresh:
                continue

            # Compute mus
            for i in range(ncenters):
                mus[i, i] = 0.0
                for j in range(i):
                    mus[i, j] = (dists[i] - dists[j]) / all_center_dists[i, j]
                    mus[j, i] = -mus[i, j]

            mus_a = mus / a
            # Eq. (14) in [1]
            # Use gg instead of g, as we already use g as index
            gg = (
                35.0 * mus_a - 35.0 * mus_a**3 + 21.0 * mus_a**5 - 5.0 * mus_a**7
            ) / 16.0
            # Apply bounds; eq. (11) in [1]
            for i in range(ncenters):
                for j in range(ncenters):
                    if i == j:
                        continue
                    mus_ij = mus[i, j]
                    if mus_ij <= -a:
                        gg[i, j] = -1.0
                    elif mus_ij >= a:
                        gg[i, j] = 1.0

            # Eq. (8) in [1]
            s_ij = 0.5 * (1.0 - gg)

            atomic_weight = np.empty(ncenters)
            for i in range(ncenters):
                atomic_weight[i] = 1.0
                for j in range(i):
                    atomic_weight[i] *= s_ij[i, j]
                for j in range(i + 1, ncenters):
                    atomic_weight[i] *= s_ij[i, j]
            part_weights[start_index + g] = atomic_weight[p] / atomic_weight.sum()


def get_mol_grid(
    atoms, coords3d, part_func=stratmann_partitioning, weight_thresh=1e-15
):
    """Get molecular grid for given atoms, centered at given Cartesian coordinates.

    Parameters
    ----------
    atoms
        List of atomic symbols.
    coords3d
        2d array denoting the origins of the atomic grids. Usually the atomic
        coordiantes.
    part_func
        Function used to repartiton the atomic grid weights.
    weight_thresh
        Positive floating point number used for grid compression. Quadrature points
        with values below this threshold will be ignored.

    Returns
    -------
    xyz
        2d array containing Cartesian gridpoints of shape (npoints, 3).
    weights
        1d array containing integration weights of shape (npoints, ).
    """
    atomic_grids = list()
    weights = list()
    for atom, origin in zip(atoms, coords3d):
        xyz, ww = get_atomic_grid(atom, origin, kind="g3")
        atomic_grids.append(xyz)
        weights.append(ww)
    part_weights = np.empty(sum([len(ag) for ag in atomic_grids]))

    numba_atomic_grids = numba.typed.List()
    for xyz in atomic_grids:
        numba_atomic_grids.append(xyz)

    part_func(coords3d, numba_atomic_grids, part_weights)

    xyz = np.concatenate(atomic_grids, axis=0)
    weights = np.concatenate(weights, axis=0)
    weights *= part_weights

    # Compress grid by dropping points below threshold.
    mask = ~(weights <= weight_thresh)
    xyz = xyz[mask]
    weights = weights[mask]
    return xyz, weights
