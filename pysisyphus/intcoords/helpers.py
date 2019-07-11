import numpy as np


def get_tangent(prims1, prims2, dihed_start, normalize=False):
    """Normalized tangent between primitive internal coordinates.

    Tangent pointing from prims1 to prims2 in primitive
    internal coordinates, taking into account the periodicity of
    dihedral angles.

    Parameters
    ----------
    prims1 : np.array
        1d-array of primitive internal coordinates in the order
        (stretches, bends, dihedrals).
    prims2 : np.array
        See prims1.
    dihed_start : int
        Index of the first primitive dihedral in prims1 and prims2.

    Returns
    -------
    tangent : np.array
        1d array containing the normalized tangent pointing from prims1 to prims2.
    """
    tangent = prims2 - prims1
    diheds = tangent[dihed_start:].copy()
    diheds_plus = diheds.copy() + 2*np.pi
    diheds_minus = diheds.copy() - 2*np.pi
    bigger = np.abs(diheds) > np.abs(diheds_plus)
    diheds[bigger] = diheds_plus[bigger]
    bigger = np.abs(diheds) > np.abs(diheds_minus)
    diheds[bigger] = diheds_minus[bigger]
    tangent[dihed_start:] = diheds
    if normalize:
        tangent /= np.linalg.norm(tangent)
    return tangent
