#!/usr/bin/env python3

def bfgs_multiply(s_list, y_list, force):
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

    r = q
    if cycles > 0:
        s = s_list[-1]
        y = y_list[-1]
        gamma = s.dot(y) / y.dot(y)
        r = gamma * q
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
        restrict_step = lambda x: x

    def lbfgs(x, *getter_args):
        nonlocal cur_cycle
        nonlocal s_list
        nonlocal y_list

        prev_forces = forces[-1]
        step = -bfgs_multiply(s_list, y_list, prev_forces)
        step = restrict_step(step)
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
