#!/usr/bin/env python3

import itertools as it

import numpy as np
import rmsd as rmsd_ext
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from pysisyphus.Geometry import Geometry


def match_geom_atoms(ref_geom, geom_to_match, hydrogen=True):
    """Apply the hungarian method to geom_to_match.

    Uses the scipy implementation of the Hungarian method/
    Kuhn-Munkres algorithm in scipy.optimize.linear_sum_assignment.
    
    See
        [1] 10.1021/ci400534h
        [2] 10.1021/acs.jcim.6b00516
    """

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
        old_inds = inds_to_match[atom]
        new_inds = old_inds[col_inds]
        new_coords_for_atom = coords_to_match_for_atom[col_inds]
        # Update coordinates and modify the array directly
        c3d = geom_matched.coords3d
        c3d[old_inds] = new_coords_for_atom
    return geom_matched


def apply_transform(coords3d, swap, reflection):
    # Swap axes
    coords3d = coords3d[:, swap]
    # Reflect
    coords3d *= reflection
    return coords3d


def match_rmsd(geom1, geom2):
    # Six possible axis swaps, (3*2, (x, y, z)*((x), (y), (z))).
    axes = (0, 1, 2)
    swaps =  list(it.permutations(axes))
    # Eight possible reflections, (±x * ±y * ±z).
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

    # Work on copies as calling standard_orientation moves the
    # geometries's coordinates.
    geom1_copy = geom1.copy()
    geom1_copy.standard_orientation()
    coords3d_1 = geom1_copy.coords3d
    geom2_copy = geom2.copy()
    geom2_copy.standard_orientation()
    coords3d_2 = geom2_copy.coords3d.copy()

    # rmsd_before = rmsd_ext.rmsd(coords3d_1, coords3d_2)
    # rmsd_before_kabsch = rmsd_ext.kabsch_rmsd(coords3d_1, coords3d_2)
    # print("Moved geometries into standard orientation.")
    # print(f"Normal RMSD before transformation: {rmsd_before:.3f}")
    # print(f"Kabsch RMSD before transformation: {rmsd_before_kabsch:.3f}")
    matched_rmsds = list()
    matched_coords = list()
    for i, transform in enumerate(transforms):
        # print(f"Cycle {i}, transformation {transform}")
        c3d = coords3d_2.copy()
        c3d_trans = apply_transform(c3d, *transform)
        # print("Applied transformation to geom_2")
        rmsd_trans = rmsd_ext.kabsch_rmsd(coords3d_1, c3d_trans)

        tmp = geom2.copy()
        tmp.coords3d = c3d_trans
        geom2_matched = match_geom_atoms(geom1_copy, tmp)
        rmsd_after = rmsd_ext.kabsch_rmsd(coords3d_1, geom2_matched.coords3d)
        matched_rmsds.append(rmsd_after)
        matched_coords.append(geom2_matched.coords)
        # print(f"RMSD transformed: {rmsd_trans:.3f}, after matching: {rmsd_after:.3f}")
        # print()
    matched_rmsds = np.array(matched_rmsds)
    min_rmsd_ind = matched_rmsds.argmin()
    min_rmsd = matched_rmsds.min()
    best_matching_coords = matched_coords[min_rmsd_ind]
    geom2_copy.coords = best_matching_coords
    # import pdb; pdb.set_trace()
    return min_rmsd, geom2_copy
