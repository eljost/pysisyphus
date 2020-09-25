import numpy as np

from pysisyphus.helpers_pure import log
from pysisyphus.intcoords import PrimTypes
from pysisyphus.intcoords import Bend


def bend_valid(coords3d, bend_ind, min_deg, max_deg):
    val = Bend._calculate(coords3d, bend_ind)
    deg = np.rad2deg(val)
    return min_deg <= deg <= max_deg


def bend_still_valid(coords3d, indices, min_deg, max_deg):
    val = Bend._calculate(coords3d, bend_ind)
    deg = np.rad2deg(val)
    # Less than, not less or equal as in "bend_valid"
    return min_deg <= deg < max_deg


def are_collinear(vec1, vec2, deg_thresh=179.5):
    thresh = np.cos(np.deg2rad(deg_thresh))
    return abs(vec1.dot(vec2)) >= abs(thresh)


def dihedral_valid(coords3d, inds, deg_thresh=179.5):
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

    valid = not (
        are_collinear(u, w, deg_thresh=deg_thresh)
        or are_collinear(v, w, deg_thresh=deg_thresh)
    )
    return valid


def dihedrals_are_valid(coords3d, dihedral_inds, logger=None):
    valid = [dihedral_valid(coords3d, inds) for inds in dihedral_inds]
    invalid = [dihedral_ind for dihedral_ind, v in zip(dihedral_inds, valid) if not v]
    if invalid:
        log(logger, f"Invalid dihedrals: {invalid}")
    all_valid = all(valid)
    return all_valid


def check_typed_prims(
    coords3d,
    typed_prims,
    bend_min_deg,
    dihed_max_deg,
    lb_min_deg,
    keep_complements=True,
):
    funcs = {
        PrimTypes.Bend: lambda indices: bend_still_valid(
            coords3d, indices, min_deg=bend_min_deg, max_deg=lb_min_deg
        ),
        PrimTypes.Torsion: lambda indices: dihedral_valid(
            coords3d, indices, deg_thresh=dihed_max_deg,
        ),
    }
    for type_, *indices in typed_prims:
        pass
    # Check bends
    # Check dihedrals
    # Check complements

    # bend_min_deg=15,
    # dihed_max_deg=175.0,
    # lb_min_deg=175.0,
