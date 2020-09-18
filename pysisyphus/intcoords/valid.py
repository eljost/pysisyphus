import numpy as np

from pysisyphus.helpers import log
from pysisyphus.intcoords import Bend


def bend_valid(coords3d, bend_ind, min_deg, max_deg):
    val = Bend._calculate(coords3d, bend_ind)
    deg = np.rad2deg(val)
    return min_deg <= deg <= max_deg


def are_collinear(vec1, vec2, thresh=1e-4):
    # ~4e-5 corresponds to 179.5°
    return 1 - abs(vec1.dot(vec2)) <= thresh


def dihedral_valid(coords3d, inds, thresh=1e-4):
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

    valid = not (are_collinear(u, w, thresh=thresh) or are_collinear(v, w, thresh=thresh))
    return valid


def dihedrals_are_valid(coords3d, dihedral_inds, logger=None):
    valid = [dihedral_valid(coords3d, inds) for inds in dihedral_inds]
    invalid = [
        dihedral_ind for dihedral_ind, v in zip(dihedral_inds, valid) if not v
    ]
    if invalid:
        log(logger, f"Invalid dihedrals: {invalid}")
    all_valid = all(valid)
    return all_valid


