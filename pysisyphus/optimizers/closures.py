# [1] https://arxiv.org/pdf/2101.04413.pdf
#     A Regularized Limited Memory BFGS method for Large-Scale Unconstrained
#     Optimization and itsefficient Implementations
#     Tankaria, Sugimoto, Yamashita, 2021
# [2] Regularization of Limited Memory Quasi-Newton Methods for Large-Scale
#     Nonconvex Minimization
#     https://arxiv.org/pdf/1911.04584.pdf
#     Kanzow, Steck 2021

from collections import deque

import numpy as np
from scipy.sparse.linalg import spsolve


def get_update_mu_reg(
    mu_min=1e-3, gamma_1=0.1, gamma_2=5.0, eta_1=0.01, eta_2=0.9, logger=None
):
    """See 5.1 in [1]"""
    assert 0.0 < mu_min
    assert 0.0 < gamma_1 <= 1.0 < gamma_2
    assert 0.0 < eta_1 < eta_2 <= 1.0

    def log(msg):
        if logger is not None:
            logger.debug(msg)

    def update_mu_reg(mu, cur_energy, trial_energy, cur_grad, step):
        """Update regularization parameter μ_reg according to ratio r.
        See 2.3 in [1].

        r = actual_reduction / predicted_reduction
        r = (f(x) - f(x + step)) / (f(x) - model_func(step, mu))

        mode_func(step, mu) = f(x) + df(x)^T step(mu) + 1/2 step(mu)^T H^-1 step(mu)

        f(x) - model_func(step, mu)
        = f(x) - (f(x) + df(x)^T step(mu) + 1/2 step(mu)^T H(mu)^-1 step(mu))
        = -df(x)^T step(mu) - 1/2 step(mu)^T H(mu)^-1 step(mu)
        = -df(x)^T step(mu) + 1/2 step(mu)^T df(x)
        = -1/2 df(x)^T step(mu)

        r = (f(x) - f(x + step) / (-1/2 df(x)^T step(mu))
        r = -2 * (f(x) - f(x + step) / (df(x)^T step(mu))

        """
        log(f"  Cur energy: {cur_energy:.8f}")
        log(f"Trial energy: {trial_energy:.8f}")
        act_change = cur_energy - trial_energy
        log(f"   Actual change: {act_change: .8f}")
        pred_change = -1 / 2 * cur_grad.dot(step)
        log(f"Predicted change: {pred_change: .8f}")
        r = act_change / pred_change
        log(f"Ratio r={r:.8f}")

        # Default case for eta_1 <= r < eta_2. Keep mu and accept step.
        mu_updated = mu
        recompute_step = False
        # Poor actual reduction. Increase shift parameter and request new step.
        if r < eta_1:
            mu_updated *= gamma_2
            recompute_step = True
            log(f"Increased μ_reg to {mu_updated:.6f}.")
        # Significant actual reduction. Reduce shift parameter und accept step.
        elif r >= eta_2:
            mu_updated = max(mu_min, mu_updated * gamma_1)
            log(f"Decreased μ_reg to {mu_updated:.6f}.")
        return mu_updated, recompute_step

    return update_mu_reg


def bfgs_multiply(
    s_list,
    y_list,
    vector,
    beta=1,
    P=None,
    logger=None,
    gamma_mult=True,
    mu_reg=None,
    inds=None,
    cur_size=None,
):
    """Matrix-vector product H·v.

    Multiplies given vector with inverse Hessian, obtained
    from repeated BFGS updates calculated from steps in 's_list'
    and gradient differences in 'y_list'.

    Based on algorithm 7.4 Nocedal, Num. Opt., p. 178."""

    assert len(s_list) == len(
        y_list
    ), "lengths of step list 's_list' and gradient list 'y_list' differ!"

    cycles = len(s_list)
    q = vector.copy()
    alphas = list()
    rhos = list()

    if mu_reg is not None:
        # Regularized L-BFGS with ŷ = y + μs, see [1]
        assert mu_reg > 0.0
        y_list_reg = list()
        for i, si in enumerate(s_list):
            yi = y_list[i]
            y_hat = yi + mu_reg * si
            if y_hat.dot(si) <= 0:
                # See 2 in [1]
                y_hat = yi + (max(0, -si.dot(yi) / si.dot(si)) + mu_reg) + si
            y_list_reg.append(y_hat)
        y_list = y_list_reg

    # Store rho and alphas as they are also needed in the second loop
    for i in reversed(range(cycles)):
        s = s_list[i]
        y = y_list[i]
        rho = 1 / y.dot(s)
        rhos.append(rho)
        try:
            alpha = rho * s.dot(q)
            q -= alpha * y
        except ValueError:
            inds_i = inds[i]
            q_ = q.reshape(cur_size, -1)
            alpha = rho * s.dot(q_[inds_i].flatten())
            # This also modifies q!
            q_[inds_i] -= alpha * y.reshape(len(inds_i), -1)
        alphas.append(alpha)

    # Restore original order, so that rho[i] = 1/s_list[i].dot(y_list[i]) etc.
    alphas = alphas[::-1]
    rhos = rhos[::-1]

    if P is not None:
        r = spsolve(P, q)
        msg = "preconditioner"
    elif gamma_mult and (cycles > 0):
        s = s_list[-1]
        y = y_list[-1]
        gamma = s.dot(y) / y.dot(y)
        r = gamma * q
        msg = f"gamma={gamma:.4f}"
    else:
        r = beta * q
        msg = f"beta={beta:.4f}"

    if mu_reg is not None:
        msg += f" and μ_reg={mu_reg:.6f}"

    if logger is not None:
        msg = f"BFGS multiply using {cycles} previous cycles with {msg}."
        if len(s_list) == 0:
            msg += "\nProduced simple SD step."
        logger.debug(msg)

    for i in range(cycles):
        s = s_list[i]
        y = y_list[i]
        try:
            beta = rhos[i] * y.dot(r)
            r += s * (alphas[i] - beta)
        except ValueError:
            inds_i = inds[i]
            r_ = r.reshape(cur_size, -1)
            beta = rhos[i] * y.dot(r_[inds_i].flatten())
            # This also modifies r!
            r_[inds_i] += s.reshape(len(inds_i), -1) * (alphas[i] - beta)

    return r


def lbfgs_closure(force_getter, M=10, beta=1, restrict_step=None):
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

        step = bfgs_multiply(s_list, y_list, force, beta=beta)
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
        dx = -beta * F

        if len(dxs) > 0:
            # Calculate gamma
            dF_F = [dF_k.dot(F) for dF_k in dFs]
            # Use least squares to avoid crashes with singular matrices
            gammas, *_ = np.linalg.lstsq(a, dF_F, rcond=None)
            gammas = gammas[:, None]
            _ = np.array(dxs) - beta * np.array(dFs)
            # Substract step correction
            dx = dx - np.sum(gammas * _, axis=0)
        dx = restrict_step(x, dx)
        dxs.append(dx)

        # Keep only informations of the last M cycles
        dxs = dxs[-M:]

        return dx, -F

    return modified_broyden


def small_lbfgs_closure(history=5, gamma_mult=True):
    """Compact LBFGS closure.

    The returned function takes two arguments: forces and prev_step.
    forces are the forces at the current iterate and prev_step is the
    previous step that lead us to the current iterate. In this way
    step restriction/line search can be done outisde of the lbfgs function.
    """

    prev_forces = None  # lgtm [py/unused-local-variable]
    grad_diffs = deque(maxlen=history)
    steps = deque(maxlen=history)
    cur_cycle = 0

    def lbfgs(forces, prev_step=None):
        nonlocal cur_cycle
        nonlocal prev_forces

        if prev_step is not None:
            steps.append(prev_step)

        # Steepest descent in the first cycle
        step = forces
        # LBFGS in the following cycles
        if cur_cycle > 0:
            grad_diffs.append(-forces - -prev_forces)
            step = bfgs_multiply(steps, grad_diffs, forces, gamma_mult=gamma_mult)

        prev_forces = forces
        cur_cycle += 1
        return step

    return lbfgs
