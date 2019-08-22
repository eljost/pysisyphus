#!/usr/bin/env python3

import autograd.numpy as np
from autograd import grad


def calc_stretch(coords3d, bond_ind):
    n, m = bond_ind
    bond = coords3d[m] - coords3d[n]
    bond_length = np.linalg.norm(bond)
    return bond_length

stretch_grad = grad(calc_stretch)


def calc_bend(coords3d, angle_ind):
    m, o, n = angle_ind
    u_dash = coords3d[m] - coords3d[o]
    v_dash = coords3d[n] - coords3d[o]
    u_norm = np.linalg.norm(u_dash)
    v_norm = np.linalg.norm(v_dash)
    u = u_dash / u_norm
    v = v_dash / v_norm
    angle_rad = np.arccos(np.dot(u, v))
    return angle_rad


bend_grad = grad(calc_bend)


def calc_dihedral(coords3d, dihedral_ind, cos_tol=1e-9):
    m, o, p, n = dihedral_ind
    u_dash = coords3d[m] - coords3d[o]
    v_dash = coords3d[n] - coords3d[p]
    w_dash = coords3d[p] - coords3d[o]
    u_norm = np.linalg.norm(u_dash)
    v_norm = np.linalg.norm(v_dash)
    w_norm = np.linalg.norm(w_dash)
    u = u_dash / u_norm
    v = v_dash / v_norm
    w = w_dash / w_norm
    phi_u = np.arccos(np.dot(u, w))
    phi_v = np.arccos(-np.dot(w, v))
    uxw = np.cross(u, w)
    vxw = np.cross(v, w)
    cos_dihed = np.dot(uxw, vxw)/(np.sin(phi_u)*np.sin(phi_v))

    # Restrict cos_dihed to [-1, 1]
    if cos_dihed >= 1 - cos_tol:
        dihedral_rad = 0
    elif cos_dihed <= -1 + cos_tol:
        dihedral_rad = np.arccos(-1)
    else:
        dihedral_rad = np.arccos(cos_dihed)

    if dihedral_rad != np.pi:
        # wxv = np.cross(w, v)
        # if wxv.dot(u) < 0:
        if vxw.dot(u) < 0:
            dihedral_rad *= -1
    return dihedral_rad


dihedral_grad = grad(calc_dihedral)
