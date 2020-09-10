import itertools as it

import numpy as np


RAD_175 = np.deg2rad(175)

class PrimInt:

    def __init__(self, inds, val, grad):
        self.inds = inds
        self.val = val
        self.grad = grad


def are_parallel(vec1, vec2, angle_ind=None, thresh=1e-6):
    dot = max(min(vec1.dot(vec2), 1), -1)
    rad = np.arccos(dot)#vec1.dot(vec2))
    # angle > 175°
    if abs(rad) > RAD_175:
        ind_str = f" ({angle_ind})" if (angle_ind is not None) else ""
    return abs(rad) > (np.pi - thresh)


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


def eval_prim_internals(cart_coords, prim_inds):
    coords3d = cart_coords.reshape(-1, 3)
    bond_inds, bend_inds, dihedral_inds = prim_inds

    def per_type(func, ind):
        val, grad = func(coords3d, ind, True)
        return PrimInt(ind, val, grad)

    bonds = [per_type(calc_stretch, ind) for ind in bond_inds]
    bends = [per_type(calc_bend, ind) for ind in bend_inds]
    dihedrals = [per_type(calc_dihedral, ind) for ind in dihedral_inds]

    return bonds, bends, dihedrals


def update_internals(new_cartesians, old_internals, prim_inds):
    prim_ints = eval_prim_internals(new_cartesians, prim_inds)
    new_internals = [prim.val for prim in it.chain(*prim_ints)]
    internal_diffs = np.array(new_internals) - old_internals

    dihedrals = prim_ints[-1]
    dihedral_num = len(dihedrals)
    dihedral_diffs = internal_diffs[-dihedral_num:]

    # Find differences that are shifted by 2*pi
    shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
    new_dihedrals = np.array([dihed.val for dihed in dihedrals])
    if any(shifted_by_2pi):
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])
        # Update values
        for dihed, new_val in zip(dihedrals, new_dihedrals):
            dihed.val = new_val

    return prim_ints


def transform_int_step(int_step, old_cart_coords, cur_internals, B_prim, prim_inds,
                       cart_rms_thresh=1e-6, logger=None):
    """Transformation is done in primitive internals, so int_step must be given
    in primitive internals and not in DLC!"""

    def log(msg):
        if logger is not None:
            logger.debug(msg)

    new_cart_coords = old_cart_coords.copy()
    remaining_int_step = int_step
    target_internals = cur_internals + int_step

    Bt_inv_prim = np.linalg.pinv(B_prim.dot(B_prim.T)).dot(B_prim)

    last_rms = 9999
    old_internals = cur_internals
    backtransform_failed = True
    for i in range(25):
        cart_step = Bt_inv_prim.T.dot(remaining_int_step)
        cart_rms = np.sqrt(np.mean(cart_step**2))
        # Update cartesian coordinates
        new_cart_coords += cart_step
        # Determine new internal coordinates
        new_prim_ints = update_internals(new_cart_coords, old_internals, prim_inds)
        new_internals = [prim.val for prim in it.chain(*new_prim_ints)]
        remaining_int_step = target_internals - new_internals
        internal_rms = np.sqrt(np.mean(remaining_int_step**2))
        log(f"Cycle {i}: rms(Δcart)={cart_rms:1.4e}, rms(Δint.) = {internal_rms:1.5e}")

        # This assumes the first cart_rms won't be > 9999 ;)
        if (cart_rms < last_rms):
            # Store results of the conversion cycle for laster use, if
            # the internal-cartesian-transformation goes bad.
            best_cycle = (new_cart_coords.copy(), new_internals.copy())
            best_cycle_ind = i
        elif i != 0:
            # If the conversion somehow fails we fallback to the best previous step.
            log(f"Backconversion failed! Falling back to step {best_cycle_ind}")
            new_cart_coords, new_internals = best_cycle
            break
        else:
            raise Exception("Internal-cartesian back-transformation already "
                            "failed in the first step. Aborting!"
            )
        old_internals = new_internals

        last_rms = cart_rms
        if cart_rms < cart_rms_thresh:
            log("Internal to cartesian transformation converged!")
            backtransform_failed = False
            break

    # if self.check_dihedrals and (not self.dihedrals_are_valid(new_cart_coords)):
        # raise NeedNewInternalsException(new_cart_coords)

    log("")

    # Return the difference between the new cartesian coordinates that yield
    # the desired internal coordinates and the old cartesian coordinates.
    cart_step = new_cart_coords - old_cart_coords
    return new_prim_ints, cart_step, backtransform_failed
