import itertools as it

import numpy as np
import rmsd as rmsd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def match_geom_atoms(ref_geom, geom_to_match, hydrogen=True):
    """Apply the hungarian method to geom_to_match.

    Uses the scipy implementation of the Hungarian method/
    Kuhn-Munkres algorithm in scipy.optimize.linear_sum_assignment.
    
    See
        [1] 10.1021/ci400534h
        [2] 10.1021/acs.jcim.6b00516
    """

    assert len(ref_geom.atoms) == len(geom_to_match.atoms), \
        "Atom numbers don't match!"

    ref_coords, _ = ref_geom.coords_by_type
    coords_to_match, inds_to_match = geom_to_match.coords_by_type
    atoms = ref_coords.keys()

    # The new geometry that will contain the perhaps reordered
    # coordinates.
    geom_matched = geom_to_match.copy()
    for atom in atoms:
        # Don't match hydrogens if not requested
        if atom.lower() == "h" and not hydrogen:
            continue
        ref_coords_for_atom = ref_coords[atom]
        coords_to_match_for_atom = coords_to_match[atom]
        # Pairwise distances between two collections
        # Atoms of ref_geom are along the rows, atoms of geom_to_match
        # along the columns.
        cd = cdist(ref_coords_for_atom, coords_to_match_for_atom)
        # Hungarian method, row_inds are returned already sorted, so
        # we only have to consider col_inds.
        row_inds, col_inds = linear_sum_assignment(cd)
        # Resort the coordinates
        new_coords_for_atom = coords_to_match_for_atom[col_inds]
        # Original indices of the coordinates in the full coords array
        # in the Geometry object.
        old_inds = inds_to_match[atom]
        # Update coordinates and modify the array directly
        c3d = geom_matched.coords3d
        c3d[old_inds] = new_coords_for_atom
    return geom_matched


def apply_transform(coords3d, swap, reflection):
    """Apply axis swaps and reflections to a coordinate array."""
    # Swap axes
    coords3d = coords3d[:, swap]
    # Reflect
    coords3d *= reflection
    return coords3d


def matched_rmsd(geom1, geom2, thresh=5e-2):
    """RMSD for optimally aligned and matched geometries.

    Returns
    -------
    matched_rmsd : float
        RMSD of optimally aligned and matched geometries.
    matched_geoms : tuple(Geometry, Geometry)
        Tuple of the optimally aligned and matched geometries.
    """

    # Work on copies of the Geometries, as calling standard_orientation
    # moves their coordinates.
    geom1_copy = geom1.copy()
    geom1_copy.standard_orientation()
    coords3d_1 = geom1_copy.coords3d
    geom2_copy = geom2.copy()
    geom2_copy.standard_orientation()
    coords3d_2 = geom2_copy.coords3d.copy()

    # After bringing the Geometries into standard orientation we may
    # still have to consider additional axis swaps and reflections to
    # allow optimal atom matching using the Hungarian method.

    # Six possible axis swaps, (3*2, (x, y, z)*((x), (y), (z))).
    axes = (0, 1, 2)
    swaps =  list(it.permutations(axes))
    # Eight possible reflections, (±x, ±y, ±z).
    reflections = (
            ( 1,  1,  1),  # no reflection
            (-1,  1,  1),  # reflect on yz plane
            ( 1, -1,  1),  # reflect on xz plane
            ( 1,  1, -1),  # reflect on xy plane
            (-1, -1,  1),  # reflect on yz and xz planes
            (-1,  1, -1),  # reflect on yz and xy planes
            ( 1, -1, -1),  # reflect on xz and xy planes
            (-1, -1, -1)   # reflect on yz, xz, xy planes
    )
    # 48 combinations of six axis swaps and eight reflections.
    transforms = list(it.product(swaps, reflections))

    matched_rmsds = list()
    matched_coords = list()
    for i, transform in enumerate(transforms):
        # Apply swap and reflection
        c3d_trans = apply_transform(coords3d_2.copy(), *transform)
        geom2_to_match = geom2.copy()
        geom2_to_match.coords3d = c3d_trans
        # Apply Hungarian method to the transformed Geometry
        geom2_matched = match_geom_atoms(geom1_copy, geom2_to_match)
        mrmsd = rmsd.kabsch_rmsd(coords3d_1, geom2_matched.coords3d)
        matched_rmsds.append(mrmsd)
        matched_coords.append(geom2_matched.coords)

        # Break when the two geometries are similar. Then we don't have to
        # apply the remaining transformations.
        if mrmsd <= thresh:
            break

    matched_rmsds = np.array(matched_rmsds)
    min_rmsd_ind = matched_rmsds.argmin()
    min_rmsd = matched_rmsds.min()
    best_matching_coords = matched_coords[min_rmsd_ind]
    geom2_copy.coords = best_matching_coords
    return min_rmsd, (geom1_copy, geom2_copy)
