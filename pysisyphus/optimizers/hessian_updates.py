#!/usr/bin/env python3

#     [1] https://link.springer.com/article/10.1007/s00214-016-1847-3
#         Birkholz, 2016

import numpy as np


def bfgs_update(H, dx, dg):
    first_term = np.outer(dg, dg) / dg.dot(dx)
    second_term = (H.dot(np.outer(dx, dx)).dot(H)
                  / dx.dot(H).dot(dx)
    )
    return first_term - second_term


def sr1_update(z, dx):
    return np.outer(z, z) / z.dot(dx)


def psb_update(z, dx):
    first_term = (np.outer(dx, z) + np.outer(z, dx)) / dx.dot(dx)
    sec_term = dx.dot(z) * np.outer(dx, dx) / dx.dot(dx)**2
    return first_term - sec_term


def flowchart_update(H, dx, dg):
    # See [1], Sec. 2, equations 1 to 3
    z = dg - H.dot(dx)
    sr1_quot = z.dot(dx) / (np.linalg.norm(z) * np.linalg.norm(dx))
    bfgs_quot = dg.dot(dx) / (np.linalg.norm(dg) * np.linalg.norm(dx))
    if sr1_quot < -0.1:
        update = sr1_update(z, dx)
        key = "SR1"
    elif bfgs_quot > 0.1:
        update = bfgs_update(H, dx, dg)
        key = "BFGS"
    else:
        update = psb_update(z, dx)
        key = "PSB"
    return update, key


def mod_flowchart_update(H, dx, dg):
    # This version seems to work too ... at least for minimizations
    # starting from a geometry near a transition state. Interesing.
    z = dg - H.dot(dx)
    quot = z.dot(dx) / (np.linalg.norm(z) * np.linalg.norm(dx))
    if quot < -0.1:
        update = bfgs_update(H, dx, dg)
        key = "BFGS"
    elif quot > 0.1:
        update = sr1_update(z, dx)
        key = "SR1"
    else:
        update = psb_update(z, dx)
        key = "PSB"
    return update, key


def bofill_update(H, dx, dg):
    z = dg - H.dot(dx)

    # Symmetric, rank-one (SR1) update
    sr1 = sr1_update(z, dx)

    # Powell (PSB) update
    powell = psb_update(z, dx)

    # Bofill mixing-factor
    mix = z.dot(dx)**2 / (z.dot(z) * dx.dot(dx))

    # Bofill update
    bofill_update = (mix * sr1) + (1 - mix)*(powell)

    return bofill_update
