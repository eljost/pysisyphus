import itertools as it

import numpy as np


class PrimInt:

    def __init__(self, inds, val, grad):
        self.inds = inds
        self.val = val
        self.grad = grad


def are_parallel(vec1, vec2, angle_ind=None, thresh=1e-6):
    dot = max(min(vec1.dot(vec2), 1), -1)
    rad = np.arccos(dot)
    return abs(rad) > (np.pi - thresh)


def calc_stretch(coords3d, bond_ind, grad=False):
    n, m = bond_ind
    bond = coords3d[m] - coords3d[n]
    bond_length = np.linalg.norm(bond)
    if grad:
        bond_normed = bond / bond_length
        row = np.zeros_like(coords3d)
        # 1 / -1 correspond to the sign factor [1] Eq. 18
        row[m,:] =  bond_normed
        row[n,:] = -bond_normed
        row = row.flatten()
        return bond_length, row
    return bond_length


def calc_bend(coords3d, angle_ind, grad=False):
    m, o, n = angle_ind
    u_dash = coords3d[m] - coords3d[o]
    v_dash = coords3d[n] - coords3d[o]
    u_norm = np.linalg.norm(u_dash)
    v_norm = np.linalg.norm(v_dash)
    u = u_dash / u_norm
    v = v_dash / v_norm
    angle_rad = np.arccos(u.dot(v))
    if grad:
        # Eq. (24) in [1]
        if are_parallel(u, v, angle_ind):
            tmp_vec = np.array((1, -1, 1))
            par = are_parallel(u, tmp_vec) and are_parallel(v, tmp_vec)
            tmp_vec = np.array((-1, 1, 1)) if par else tmp_vec
            w_dash = np.cross(u, tmp_vec)
        else:
            w_dash = np.cross(u, v)
        w_norm = np.linalg.norm(w_dash)
        w = w_dash / w_norm
        uxw = np.cross(u, w)
        wxv = np.cross(w, v)

        row = np.zeros_like(coords3d)
        #                  |  m  |  n  |  o  |
        # -----------------------------------
        # sign_factor(amo) |  1  |  0  | -1  | first_term
        # sign_factor(ano) |  0  |  1  | -1  | second_term
        first_term = uxw / u_norm
        second_term = wxv / v_norm
        row[m,:] = first_term
        row[o,:] = -first_term - second_term
        row[n,:] = second_term
        row = row.flatten()
        return angle_rad, row
    return angle_rad


def calc_dihedral(coords3d, dihedral_ind, grad=False, cos_tol=1e-9):
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
    phi_u = np.arccos(u.dot(w))
    phi_v = np.arccos(-w.dot(v))
    uxw = np.cross(u, w)
    vxw = np.cross(v, w)
    cos_dihed = uxw.dot(vxw)/(np.sin(phi_u)*np.sin(phi_v))

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
    if grad:
        row = np.zeros_like(coords3d)
        #                  |  m  |  n  |  o  |  p  |
        # ------------------------------------------
        # sign_factor(amo) |  1  |  0  | -1  |  0  | 1st term
        # sign_factor(apn) |  0  | -1  |  0  |  1  | 2nd term
        # sign_factor(aop) |  0  |  0  |  1  | -1  | 3rd term
        # sign_factor(apo) |  0  |  0  | -1  |  1  | 4th term
        sin2_u = np.sin(phi_u)**2
        sin2_v = np.sin(phi_v)**2
        first_term  = uxw/(u_norm*sin2_u)
        second_term = vxw/(v_norm*sin2_v)
        third_term  = uxw*np.cos(phi_u)/(w_norm*sin2_u)
        fourth_term = -vxw*np.cos(phi_v)/(w_norm*sin2_v)
        row[m,:] = first_term
        row[n,:] = -second_term
        row[o,:] = -first_term + third_term - fourth_term
        row[p,:] = second_term - third_term + fourth_term
        row = row.flatten()
        return dihedral_rad, row
    return dihedral_rad


def eval_prim_internals(coords3d, prim_inds):
    bond_inds, bend_inds, dihedral_inds = prim_inds

    def per_type(func, ind):
        val, grad = func(coords3d, ind, True)
        return PrimInt(ind, val, grad)

    bonds = [per_type(calc_stretch, ind) for ind in bond_inds]
    bends = [per_type(calc_bend, ind) for ind in bend_inds]
    dihedrals = [per_type(calc_dihedral, ind) for ind in dihedral_inds]

    return bonds, bends, dihedrals


def eval_B(coords3d, prim_inds):
    prim_internals = eval_prim_internals(coords3d, prim_inds)
    return np.array([prim.grad for prim in it.chain(*prim_internals)])
