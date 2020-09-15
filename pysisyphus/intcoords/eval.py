import logging
import random

import numpy as np


class PrimInternal:

    def __init__(self, inds, val, grad):
        self.inds = inds
        self.val = val
        self.grad = grad


# PrimInternal = namedtuple("PrimitiveInternal", "inds val grad")


def are_parallel(u, v, thresh=1e-6):
    dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return (1 - abs(dot)) < thresh


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


# def eval_prim_internals(coords3d, primitives):
    # bond_inds, bend_inds, dihedral_inds = prim_inds

    # def per_type(func, ind):
        # val, grad = func(coords3d, ind, True)
        # return PrimInt(ind, val, grad)

    # bonds = [per_type(calc_stretch, ind) for ind in bond_inds]
    # bends = [per_type(calc_bend, ind) for ind in bend_inds]
    # dihedrals = [per_type(calc_dihedral, ind) for ind in dihedral_inds]

    # return bonds, bends, dihedrals


def eval_primitives(coords3d, primitives):
    prim_internals = list()
    for primitive in primitives:
        value, gradient = primitive.calculate(coords3d, gradient=True)
        prim_internal = PrimInternal(primitive.indices, value, gradient)
        prim_internals.append(prim_internal)
    return prim_internals


def eval_B(coords3d, primitives):
    prim_internals = eval_primitives(coords3d, primitives)
    return np.array([prim_int.grad for prim_int in prim_internals])


def check_primitives(coords3d, primitives, thresh=1e-6, logger=None):
    def log(msg, level=logging.DEBUG):
        if logger is not None:
            logger.log(level, msg)

    B = eval_B(coords3d, primitives)
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
