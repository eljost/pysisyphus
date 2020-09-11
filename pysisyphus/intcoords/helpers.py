import logging
import random

import numpy as np

from pysisyphus.intcoords import RedundantCoordsV2 as RedundantCoords
from pysisyphus.intcoords.eval import eval_prim_internals, eval_B


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


def get_ind_sets(geom, strict=False):
    """Convert RedundandCoords.prim_indices to sets of tuples."""
    try:
        bonds, bends, dihedrals = geom.internal.prim_indices
    # Exception will be raised if the geom is in cartesians. Then we
    # just create the internals for the given geometry.
    except AttributeError as err:
        if strict:
            raise err
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
    def not_cartesian(geom):
        return geom.coord_type != "cart"

    assert not_cartesian(geom1) and not_cartesian(geom2)
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


def are_collinear(vec1, vec2, thresh=1e-4):
    # ~4e-5 corresponds to 179.5°
    return 1 - abs(vec1.dot(vec2)) <= thresh


def dihedrals_are_valid(cart_coords, dihedral_inds):
    coords3d = cart_coords.reshape(-1, 3)

    def valid_dihedral(inds):
        m, o, p, n = inds
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[p]
        w_dash = coords3d[p] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w_norm = np.linalg.norm(w_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        w = w_dash / w_norm

        valid = not (are_collinear(u, w) or are_collinear(v, w))
        return valid

    all_valid = all([valid_dihedral(inds) for inds in dihedral_inds])
    return all_valid


def check_primitives(coords3d, prim_inds, thresh=1e-6, logger=None):
    def log(msg, level=logging.DEBUG):
        if logger is not None:
            logger.log(level, msg)

    B = eval_B(coords3d, prim_inds)
    G = B.T.dot(B)
    w, v = np.linalg.eigh(G)
    nonzero_inds = np.abs(w) > thresh
    # More coordinates may be expected when collinear atoms are present.
    expected = coords3d.size - 6
    nonzero_num = sum(nonzero_inds)
    missing = expected - nonzero_num
    if missing > 0:
        log( "Not enough internal coordinates defined! Expected at least "
            f"{expected} nonzero eigenvalues, but found only {nonzero_num}!"
        )
    nonzero_w = w[nonzero_inds]
    # Condition number
    kappa = abs(nonzero_w.max()/nonzero_w.min())
    log(f"Condition number of B^T.B=G: {kappa:.2e}")
    return missing+1, kappa


def augment_primitives(missing_prims, coords3d, prim_indices, fragments):
    add_bonds = list()
    add_bends = list()
    add_dihedrals = list()

    fragment_tpls = [tuple(fragment) for fragment in fragments]
    if len(fragments) > 1:
        frag_inds = list(range(len(fragments)))
        bond_inds = prim_indices[0]
        bond_sets = [frozenset(bond) for bond in bond_inds]
        while missing_prims > 0:
            random.shuffle(fragment_tpls)
            frag1, frag2 = fragment_tpls[:2]
            atom1 = random.choice(frag1)
            atom2 = random.choice(frag2)
            bond_set = frozenset((atom1, atom2))
            if (bond_set not in bond_sets) and (bond_set not in add_bonds):
                add_bonds.append(list(bond_set))
                missing_prims -= 1
    return add_bonds, add_bends, add_dihedrals
