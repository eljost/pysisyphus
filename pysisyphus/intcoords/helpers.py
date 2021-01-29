import itertools as it

import numpy as np

from pysisyphus.intcoords.PrimTypes import PrimTypes, PrimTypeShortcuts
from pysisyphus.intcoords import RedundantCoords
from pysisyphus.intcoords.exceptions import DifferentPrimitivesException


def get_tangent(prims1, prims2, dihedral_inds, normalize=False):
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
    dihedral_inds : list of int
        Dihedral indices in prims1 and prims2.

    Returns
    -------
    tangent : np.array
        1d array containing the normalized tangent pointing from prims1 to prims2.
    """
    try:
        tangent = prims2 - prims1
    except ValueError:
        raise DifferentPrimitivesException
    diheds = tangent[dihedral_inds].copy()
    diheds_plus = diheds.copy() + 2 * np.pi
    diheds_minus = diheds.copy() - 2 * np.pi
    bigger = np.abs(diheds) > np.abs(diheds_plus)
    diheds[bigger] = diheds_plus[bigger]
    bigger = np.abs(diheds) > np.abs(diheds_minus)
    diheds[bigger] = diheds_minus[bigger]
    tangent[dihedral_inds] = diheds
    if normalize:
        tangent /= np.linalg.norm(tangent)
    return tangent


def get_step(geom, coords):
    assert len(geom.coords) == len(coords)

    if geom.coord_type == "cart":
        diff = geom.coords - coords
    elif geom.coord_type in ("redund", "dlc"):
        diff = -get_tangent(
            geom.internal.prim_coords, coords, geom.internal.dihedral_inds
        )
    else:
        raise Exception("Invalid coord_type!")

    # Convert to DLC
    if geom.coord_type == "dlc":
        diff = geom.internal.U.T.dot(diff)

    return diff


def merge_coordinate_definitions(geom1, geom2):
    typed_prims1 = geom1.internal.typed_prims
    typed_prims2 = geom2.internal.typed_prims
    union = list(set(typed_prims1) | set(typed_prims2))
    union.sort()
    # Check if internal coordinates that are only present in one of the two
    # geometries are valid in the other. If not, we omit these primitives.
    redundant = RedundantCoords(geom1.atoms, geom1.cart_coords, typed_prims=union)
    valid_typed_prims = redundant.typed_prims
    return valid_typed_prims


def form_coordinate_union(geom1, geom2):
    def not_cartesian(geom):
        return geom.coord_type != "cart"

    assert not_cartesian(geom1) and not_cartesian(geom2)
    # The first call yields all primitives from geom1 that are also valid at geom2.
    typed_prims1 = merge_coordinate_definitions(geom1, geom2)
    # The second call yields all primitives from geom2 that are also valid at geom1.
    typed_prims2 = merge_coordinate_definitions(geom2, geom1)

    intersection = list(set(typed_prims1) & set(typed_prims2))
    intersection.sort()
    return intersection


def normalize_prim_input(prim_inp):
    """Normalize input for define_prims and constrain_prims

    The intcoords.RedundantCoords constructor expects lists of integer lists
    (tuples) for arguments like 'define_prims' and 'constrain_prims'. The first item
    of every list determines the type of primitive coordinate. Currently
    there are about 20 different types and it is hard to remember all of
    them.

    So we also allow a more human friendly input, that is normalized here.
    The most common primitives are:

    0: BOND
    5: BEND
    8: PROPER_DIHEDRAL

    This function maps inputs like ["BOND", 1, 2] to [PrimTypes.BOND, 1, 2] etc.

    Always returns a list of tuples, as some prim_inps expand to multiple
    coordinates, e.g., XYZ or ATOM.
    """
    prim_type, *indices = prim_inp

    # First check if we got something like an integer
    try:
        return [tuple([PrimTypes(int(prim_type))] + indices)]
    # Raised when prim_type is, e.g., "BOND"
    except ValueError:
        pass

    # Check if we got a PrimType name
    try:
        prim_type_ = getattr(PrimTypes, str(prim_type).upper())
        return [tuple([prim_type_] + indices)]
    except AttributeError:
        pass

    # Check if we got a shortcut, e.g, X/Y/Z/XYZ/ATOM etc.
    try:
        prim_types_ = PrimTypeShortcuts[str(prim_type).upper()]
        return [tuple([prim_type_] + indices) for prim_type_ in prim_types_]
    except KeyError as error:
        print(f"Could not normalize 'prim_inp'={prim_inp}!")
        raise error


def normalize_prim_inputs(prim_inps):
    # Flatten list of tuples
    return list(it.chain(*[normalize_prim_input(pi) for pi in prim_inps]))
