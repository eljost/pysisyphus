#!/usr/bin/env python3

# [1] https://link.springer.com/article/10.1007/s00214-016-1847-3
#     Birkholz, 2016
# [2] Geometry optimization in Cartesian coordinates: Constrained optimization
#     Baker, 1992
# [3] https://epubs.siam.org/doi/pdf/10.1137/S1052623496306450
#     BFGS WITH UPDATE SKIPPING AND VARYING MEMORY
#     Kolda, 1998
# [4] https://link.springer.com/article/10.1186/1029-242X-2012-241
#     New cautious BFGS algorithm based on modified Armijo-type line search
#     Wan, 2012
# [5] Numerical optimization, 2nd ed.
#     Nocedal, Wright

import numpy as np


def bfgs_update(H, dx, dg):
    first_term = np.outer(dg, dg) / dg.dot(dx)
    second_term = (H.dot(np.outer(dx, dx)).dot(H)
                  / dx.dot(H).dot(dx)
    )
    return first_term - second_term, "BFGS"


def damped_bfgs_update(H, dx, dg):
    """See [5]"""
    dxdg = dx.dot(dg)
    dxHdx = dx.dot(H).dot(dx)
    theta = 1
    if dxdg < 0.2*dxHdx:
        theta = 0.8*dxHdx / (dxHdx - dxdg)
    r = theta*dg + (1-theta)*H.dot(dx)

    first_term = np.outer(r, r) / r.dot(dx)
    second_term = (H.dot(np.outer(dx, dx)).dot(H)
                  / dxHdx
    )
    return first_term - second_term, "damped BFGS"


def sr1_update(z, dx):
    return np.outer(z, z) / z.dot(dx), "SR1"


def psb_update(z, dx):
    first_term = (np.outer(dx, z) + np.outer(z, dx)) / dx.dot(dx)
    sec_term = dx.dot(z) * np.outer(dx, dx) / dx.dot(dx)**2
    return first_term - sec_term, "PSB"


def flowchart_update(H, dx, dg):
    # See [1], Sec. 2, equations 1 to 3
    z = dg - H.dot(dx)
    sr1_quot = z.dot(dx) / (np.linalg.norm(z) * np.linalg.norm(dx))
    bfgs_quot = dg.dot(dx) / (np.linalg.norm(dg) * np.linalg.norm(dx))
    if sr1_quot < -0.1:
        update, key = sr1_update(z, dx)
    elif bfgs_quot > 0.1:
        update, key = bfgs_update(H, dx, dg)
    else:
        update, key = psb_update(z, dx)
    return update, key


def mod_flowchart_update(H, dx, dg):
    # This version seems to work too ... at least for minimizations
    # starting from a geometry near a transition state. Interesing.
    z = dg - H.dot(dx)
    quot = z.dot(dx) / (np.linalg.norm(z) * np.linalg.norm(dx))
    if quot < -0.1:
        update, key = bfgs_update(H, dx, dg)
    elif quot > 0.1:
        update, key = sr1_update(z, dx)
    else:
        update, key = psb_update(z, dx)
    return update, key


def bofill_update(H, dx, dg):
    z = dg - H.dot(dx)

    # Symmetric, rank-one (SR1) update
    sr1, _ = sr1_update(z, dx)

    # Powell (PSB) update
    powell, _ = psb_update(z, dx)

    # Bofill mixing-factor
    mix = z.dot(dx)**2 / (z.dot(z) * dx.dot(dx))

    # Bofill update
    bofill_update = (mix * sr1) + (1 - mix)*(powell)

    return bofill_update, "Bofill"


"""
def multi_step_update(H, coords, gradients, energies, last_cycles=3,
                      key="flowchart"):
    coords = np.array(coords)
    gradients = np.array(gradients)
    energies = np.array(energies)

    # At the beginning there may not be enough cycles present.
    # Substract 1 because N coords are produced from N-1 steps.
    use_cycles = min(len(coords)-1, last_cycles)
    funcs = {
        "bfgs": bfgs_update,
        "flowchart": flowchart_update,
        "bofill": bofill_update, 
    }
    update_func = funcs[key]
    steps = np.diff(coords, axis=0)[-last_cycles:]
    grad_diffs = np.diff(gradients, axis=0)[-last_cycles:]
    energy_diffs = np.diff(energies)[-last_cycles:]
    # Exclude steps that yielded and increased energy
    energy_lowered = energy_diffs < 0
    # See [2] p. 247, step 9, regarding the hessian update.
    # Exclude cycles that would lead to positive but small divisors
    # and exclude cycles where the dot product is negative (gradient
    # and step vector point in opposite directions, 'curvature conditon').
    # The einsum call calcualtes the row-wise dot product
    dgdx_big_enough = np.einsum("ij,ij->i", grad_diffs, steps) > 1e-8
    use_cycles = energy_lowered & dgdx_big_enough
    last_coord = coords[0]
    last_grad = gradients[0]
    for coord, grad in zip(coords, gradients):

    # Recalculate dx and dg values because we maybe left out some cycles
    dxs = list()
    dgs = list()
    for i, use in enumerate(use_cycles):
        if use:
            dxs.append(steps[-last_cycles+i])
            dgs.append(grad_diffs[-last_cycles+i])
    import pdb; pdb.set_trace()
    grads_ = gradients[-last_cycles:][use_cycles]
    coords_ = coords[-last_cycles:][use_cycles]
    grad_diffs_ = np.diff(grads_, axis=0)
    steps_ = np.diff(coords_, axis=0)
    for dx, dg in zip(steps_, grad_diffs_):
        dH = update_func(H, dx, dg)
        H += dH
    return H
"""


# def multi_step_update(H, coords, gradients, energies, last_cycles=3,
                      # key="bfgs", ref=None, steps_ref=None):
def multi_step_update(H, steps, gradients, energies, last_cycles=1, key="bfgs"):
    last_cycles = min(len(steps), last_cycles)
    funcs = {
        "bfgs": bfgs_update,
        "flowchart": flowchart_update,
        "bofill": bofill_update, 
    }
    update_func = funcs[key]
    # TODO: calculation of steps is wrong/not so easy ;)
    steps_from_current = np.cumsum(steps, axis=0)
    grad_diffs = gradients[-1] - np.array(gradients)[-last_cycles-1:-1]
    energy_diffs = energies[-1] - np.array(energies)[-last_cycles-1:-1]
    # Exclude steps that yielded and increased energy
    energy_lowered = energy_diffs < 0
    # See [2] p. 247, step 9, regarding the hessian update.
    # Exclude cycles that would lead to positive but small divisors
    # and exclude cycles where the dot product is negative (gradient
    # and step vector point in opposite directions, 'curvature conditon').
    # The einsum call calcualtes the row-wise dot product
    # dgdx_big_enough = np.einsum("ij,ij->i", grad_diffs, steps) > 1e-8
    # use_cycles = energy_lowered & dgdx_big_enough
    # use_cycles = energy_lowered
    dgdxs = np.einsum("ij,ij->i", grad_diffs, steps_from_current)
    # dgdx_big_enough = np.einsum("ij,ij->i", grad_diffs, steps) > 1e-8
    dgdx_big_enough = np.abs(dgdxs) > 1e-8
    use_cycles = dgdx_big_enough
    # print("Use cycles:", use_cycles)
    H_updated = H.copy()
    if len(gradients) == 3:
        import pdb; pdb.set_trace()
    for i, to_use in enumerate(use_cycles[::-1], 1):
        if not to_use:
            continue
        dx = steps_from_current[-i]
        dg = grad_diffs[-i]
        dH, _ = update_func(H, dx, dg)
        H_updated += dH
    return H_updated
