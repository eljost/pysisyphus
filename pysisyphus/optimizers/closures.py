#!/usr/bin/env python3

import numpy as np


def bfgs_multiply(s_list, y_list, force, beta=1):
    """Get a L-BFGS step.
    
    Algorithm 7.4 Nocedal, Num. Opt., p. 178."""
    q = -force
    cycles = len(s_list)
    alphas = list()
    rhos = list()
    # Store rho and alphas as they are also needed in the second loop
    for i in reversed(range(cycles)):
        s = s_list[i]
        y = y_list[i]
        rho = 1/y.dot(s)
        rhos.append(rho)
        alpha = rho * s.dot(q)
        alphas.append(alpha)
        q = q - alpha*y
    # Restore original order, so that rho[i] = 1/s_list[i].dot(y_list[i]) etc.
    alphas = alphas[::-1]
    rhos = rhos[::-1]

    if cycles > 0:
        s = s_list[-1]
        y = y_list[-1]
        gamma = s.dot(y) / y.dot(y)
        r = gamma * q
    else:
        r = beta * q
    for i in range(cycles):
        s = s_list[i]
        y = y_list[i]
        beta = rhos[i] * y.dot(r)
        r = r + s*(alphas[i] - beta)

    return r


def lbfgs_closure(first_force, force_getter, m=10, restrict_step=None):
    s_list = list()
    y_list = list()
    forces = [first_force, ]
    cur_cycle = 0

    if restrict_step is None:
        restrict_step = lambda x, dx: dx

    def lbfgs(x, *getter_args):
        nonlocal cur_cycle
        nonlocal s_list
        nonlocal y_list

        prev_forces = forces[-1]
        step = -bfgs_multiply(s_list, y_list, prev_forces)
        step = restrict_step(x, step)
        new_x = x + step
        new_forces = force_getter(new_x, *getter_args)
        s = new_x - x
        s_list.append(s)
        y = prev_forces - new_forces
        y_list.append(y)
        forces.append(new_forces)
        # Only keep last m cycles
        s_list = s_list[-m:]
        y_list = y_list[-m:]
        cur_cycle += 1
        return new_x, step, new_forces
    return lbfgs


def lbfgs_closure_(force_getter, M=10, beta=1, restrict_step=None):
    x_list = list()
    s_list = list()
    y_list = list()
    force_list = list()
    cur_cycle = 0

    if restrict_step is None:
        restrict_step = lambda x, dx: dx

    def lbfgs(x, *getter_args):
        nonlocal x_list
        nonlocal cur_cycle
        nonlocal s_list
        nonlocal y_list

        force = force_getter(x, *getter_args)
        if cur_cycle > 0:
            prev_x = x_list[-1]
            s = x - prev_x
            s_list.append(s)
            prev_force = force_list[-1]
            y = prev_force - force
            y_list.append(y)
        x_list.append(x)
        force_list.append(force)

        step = -bfgs_multiply(s_list, y_list, force, beta=beta)
        step = restrict_step(x, step)
        # Only keep last m cycles
        s_list = s_list[-M:]
        y_list = y_list[-M:]
        cur_cycle += 1
        return step, force
    return lbfgs


def modified_broyden_closure(force_getter, M=5, beta=1, restrict_step=None):
    """https://doi.org/10.1006/jcph.1996.0059
    F corresponds to the residual gradient, so we after calling
    force_getter we multiply the force by -1 to get the gradient."""

    dxs = list()
    dFs = list()
    # As used in scipy
    x = None
    F = None
    # The next line is used in SciPy, but then beta strongly depends
    # on the magnitude of x, which is quite weird. Disregarding the
    # analytical potentials the electronic energy is usually
    # invariant under translation and rotation and doesn't depend on
    # the magnitude of x.
    # beta = 0.5 * max(np.linalg.norm(x), 1) / np.linalg.norm(F)
    a = None

    if restrict_step is None:
        restrict_step = lambda x, dx: dx

    def modified_broyden(x, *getter_args):
        nonlocal dxs
        nonlocal dFs
        nonlocal F
        nonlocal beta
        nonlocal a

        F_new = -force_getter(x, *getter_args)
        if F is not None:
            dF = F_new - F
            dFs.append(dF)
            dFs = dFs[-M:]
            # Overlap matrix
            a = np.zeros((len(dFs), len(dFs)))
            for k, dF_k in enumerate(dFs):
                for m, dF_m in enumerate(dFs):
                    a[k, m] = dF_k.dot(dF_m)
        F = F_new
        dx = -beta*F

        if len(dxs) > 0:
            # Calculate gamma
            dF_F = [dF_k.dot(F) for dF_k in dFs]
            gammas = np.linalg.solve(a, dF_F)[:,None]
            _ = np.array(dxs) - beta*np.array(dFs)
            # Substract step correction
            dx = dx - np.sum(gammas * _, axis=0)
        dx = restrict_step(x, dx)
        x_new = x + dx
        dxs.append(dx)

        # Keep only informations of the last M cycles
        dxs = dxs[-M:]

        return dx, -F
    return modified_broyden
