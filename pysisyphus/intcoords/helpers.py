#!/usr/bin/env python3

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords


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


def get_step(geom, coords):
    assert len(geom.coords) == len(coords)

    if geom.coord_type == "cart":
        diff = geom.coords - coords
    elif geom.coord_type in ("redund", "dlc"):
        diff = -get_tangent(geom.internal.prim_coords, coords,
                            geom.internal.dihed_start)
    else:
        raise Exception("Invalid coord_type!")

    # Convert to DLC
    if geom.coord_type == "dlc":
        diff = geom.internal.U.T.dot(diff)

    return diff


def to_set(iterable):
    """Convert iterable of iterable to a set of tuples."""
    return set([tuple(_) for _ in iterable])


def get_ind_sets(geom):
    """Convert RedundandCoords.prim_indices to sets of tuples."""
    try:
        bonds, bends, dihedrals = geom.internal.prim_indices
    # Exception will be raised if the geom is in cartesians. Then we
    # just create the internals for the given geometry.
    except AttributeError:
        internal = RedundantCoords(geom.atoms, geom.cart_coords)
        bonds, bends, dihedrals = internal.prim_indices
    return to_set(bonds), to_set(bends), to_set(dihedrals)


def merge_coordinate_definitions(geom1, geom2):
    bonds1, bends1, dihedrals1 = get_ind_sets(geom1)
    bonds2, bends2, dihedrals2 = get_ind_sets(geom2)
    # Form new superset of coordinate definitions that contain
    # all definitions from geom1 and geom2.
    all_bonds = bonds1 | bonds2
    all_bends = bends1 | bends2
    all_dihedrals = dihedrals1 | dihedrals2
    all_prim_indices = (all_bonds, all_bends, all_dihedrals)
    # Check if internal coordinates that are only present in one
    # of the two geometries are valid in the other. If not we omit
    # this coordinate definition in the end.
    redundant = RedundantCoords(geom1.atoms, geom1.cart_coords,
                                prim_indices=all_prim_indices)
    bonds, bends, dihedrals = redundant.prim_indices
    return to_set(bonds), to_set(bends), to_set(dihedrals)


def form_coordinate_union(geom1, geom2):
    # The first call yields all primitives from geom1 that are also
    # valid at geom2.
    bonds1, bends1, dihedrals1 = merge_coordinate_definitions(geom1, geom2)
    # The second call yields all primitives from geom2 that are also
    # valid at geom1.
    bonds2, bends2, dihedrals2 = merge_coordinate_definitions(geom2, geom1)

    # Only use primitive coordinate definitions that are valid for both
    valid_bonds = bonds1 & bonds2
    valid_bends = bends1 & bends2
    valid_dihedrals = dihedrals1 & dihedrals2
    prim_indices = (valid_bonds, valid_bends, valid_dihedrals)
    # return valid_bonds, valid_bends, valid_dihedrals
    return prim_indices
