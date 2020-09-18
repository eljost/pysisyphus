# [1] Nocedal - Numerical Optimization, Second ed.

import logging

import numpy as np


logger = logging.getLogger("optimizer")
def log(msg):
    logger.debug(msg)


def interpol_alpha_quad(f_0, df_0, f_alpha_0, alpha_0):
    return -df_0*alpha_0**2 / 2 / (f_alpha_0 - f_0 - df_0*alpha_0)


def interpol_alpha_cubic(f_0, df_0, f_alpha_0, f_alpha_1, alpha_0, alpha_1):
    quot = 1 / (alpha_0**2 * alpha_1**2 * (alpha_1 - alpha_0))
    A = np.array(((alpha_0**2, -alpha_1**2),
                  (-alpha_0**3, alpha_1**3)))
    B = np.array(((f_alpha_1 - f_0 - df_0*alpha_1),
                  (f_alpha_0 - f_0 - df_0*alpha_0)))
    a, b = quot * A @ B
    alpha_cubic = (-b + (b**2 - 3*a*df_0)**(0.5)) / (3*a)
    return alpha_cubic


class LineSearchConverged(Exception):

    def __init__(self, alpha):
        self.alpha = alpha


class LineSearchNotConverged(Exception):
    pass


def linesearch_wrapper(cond):
    def cond_wrapper(func):
        def wrapper(f, df, x0, p, f0=None, g0=None, c1=0.1, c2=0.9, max_cycles=10,
                    *args, **kwargs):
            alpha_fs = {}
            def _phi(alpha):
                alpha = float(alpha)
                try:
                    f_alpha = alpha_fs[alpha]
                except KeyError:
                    log(f"\tEvaluating energy for alpha={alpha:.6f}")
                    f_alpha = f(x0 + alpha*p)
                    alpha_fs[alpha] = f_alpha
                return f_alpha

            alpha_gs = {}
            dphis = {}
            def _dphi(alpha):
                alpha = float(alpha)
                try:
                    df_alpha = alpha_gs[alpha]
                    dphi_ = df_alpha @ p
                except KeyError:
                    log(f"\tEvaluating gradient for alpha={alpha:.6f}")
                    df_alpha = df(x0 + alpha*p)
                    alpha_gs[alpha] = df_alpha
                    dphi_ = df_alpha @ p
                    dphis[alpha] = dphi_
                return dphi_

            def get_phi_dphi(what, alpha, check=True):
                """Wrapper that handles function/gradient evaluations."""
                alpha = float(alpha)
                whats = "f g fg".split()
                assert what in whats
                calc_funcs = {
                    "f": _phi,
                    "g": _dphi,
                }
                result = [calc_funcs[w](alpha) for w in what]
                # Check if we got both phi and dphi for alpha now. If so we
                # can check if the chosen condition (Wolfe/approx. Wolfe) is
                # satisfied.
                if check and (alpha > 0.0) \
                   and (alpha in alpha_fs) and (alpha in alpha_gs) and cond_func(alpha):
                    raise LineSearchConverged(alpha)
                # Dont return a list if only f or g was requested.
                if len(what) == 1:
                    result = result[0]
                return result

            def get_fg(what, alpha):
                """Lookup raw function/gradient values at alpha."""
                whats = "f g fg".split()
                assert what in whats
                lookups = {
                    "f": alpha_fs,
                    "g": alpha_gs,
                }
                result = [lookups[w][alpha] for w in what]
                if len(what) == 1:
                    result = result[0]
                return result

            if f0 is None:
                phi0 = get_phi_dphi("f", 0)
            else:
                phi0 = f0
                alpha_fs[0.] = f0
            if g0 is None:
                dphi0 = get_phi_dphi("g", 0)
            else:
                dphi0 = g0 @ p
                alpha_gs[0.] = g0

            def sufficiently_decreased(alpha):
                """Sufficient decrease/Armijo condition."""
                return _phi(alpha) <= (phi0 + c1 * alpha * dphi0)

            def curvature_condition(alpha):
                return _dphi(alpha) >= c2*dphi0

            def strong_curvature_condition(alpha):
                return abs(_dphi(alpha)) <= -c2*dphi0

            def wolfe_condition(alpha):
                """Normal, not strong, Wolfe condition."""
                return sufficiently_decreased(alpha) \
                       and curvature_condition(alpha)

            def strong_wolfe_condition(alpha):
                """Strong wolfe condition."""
                return sufficiently_decreased(alpha) \
                       and strong_curvature_condition(alpha)

            conds = {
                "armijo": sufficiently_decreased,
                "curv": curvature_condition,
                "wolfe": wolfe_condition,
                "strong_wolfe": strong_wolfe_condition,
            }

            cond_func = conds[cond]

            # Q_prev = 0
            # C_prev = 0
            # print(f"\talpha_init={alpha_init:.6f}, ak={ak:.6f}, bk={bk:.6f}")
            # approx_permanently = False
            # conds = {
                # # True: ("approx.", t2_condition),
                # True: ("approx.", approx_wolfe_condition),
                # False: ("wolfe", wolfe_condition),
            # }
            # if f_prev:
                # Q = 1 + Q_prev * Delta
                # C = C_prev + (abs(alpha_fs[0.]) - C_prev)/Q
                # use_approx = abs(f0 - f_prev) <= omega*C
                # if approx_permanently or use_approx:
                    # approx_permanently = True
            # else:
                # use_approx = False
            # cond_name, cond = conds[use_approx]
            # # print(f"{k:02d}: [{ak:.6f},{bk:.6f}]")
            # # if wolfe_condition(ak) or t2_condition(ak):
            # print(f"Using {cond_name} condition")
            # cond = wolfe_condition
            # cond = t2_condition
            # cond = approx_wolfe_condition

            linesearch_result = func(x0, p, get_phi_dphi, get_fg, conds,
                                     max_cycles, *args, **kwargs)
            return linesearch_result
        return wrapper
    return cond_wrapper


@linesearch_wrapper("wolfe")
def hager_zhang(x0, p, get_phi_dphi, get_fg, conds, max_cycles,
                alpha_init=None, alpha_prev=None,
                f_prev=None, dphi0_prev=None, quad_step=False,
                eps=1e-6, theta=0.5, gamma=0.5, rho=5,
                psi_0=.01, psi_1=.1, psi_2=2., psi_low=0.1, psi_hi=10,
                Delta=.7, omega=1e-3, max_bisects=10):
    epsk = eps * abs(get_fg("f", 0.))
    phi0, dphi0 = get_phi_dphi("fg", 0.)
    f0, g0 = get_fg("fg", 0.)

    cond = conds["wolfe"]

    import pdb; pdb.set_trace()
    def bisect(a, b):
        """Bisect interval [a, b]."""
        for i in range(max_bisects):
            # U3 a.
            d = (1 - theta)*a + theta*b
            dphi_d = get_phi_dphi("g", d)
            if dphi_d >= 0:
                return a, d

            phi_d = get_phi_dphi("f", d)
            # U3 b.
            # If (dphi_d > 0) we would already have returned above...
            if phi_d <= phi0 + epsk:
                a = d
            # U3 c.
            elif phi_d > phi0 + epsk:
                b = d
        raise Exception("Bisect failed!")

    def interval_update(a, b, c):
        """Narrows down the bracketing interval."""
        # U0
        if not (a < c < b):
            return a, b

        phi_c, dphi_c = get_phi_dphi("fg", c)
        # U1, sign of slope projection changed. We already passed the minimum.
        if dphi_c >= 0:
            return a, c
        # U2, we are moving towards the minimum.
        elif phi_c <= phi0 + epsk:
            return c, b

        # U3, phi_c increased above phi0, so we probably passed the minimum.
        return bisect(a, c)

    def secant(a, b):
        """Take secant step."""
        dphia = get_phi_dphi("g", a)
        dphib = get_phi_dphi("g", b)
        return (a*dphib - b*dphia) / (dphib - dphia)

    def double_secant(a, b):
        """Take secant² step."""
        c = secant(a, b)
        A, B = interval_update(a, b, c)
        cB_close = np.isclose(c, B)
        cA_close = np.isclose(c, A)

        if cB_close:
            c_dash = secant(b, B)
        elif cA_close:
            c_dash = secant(a, A)

        if cB_close or cA_close:
            a_dash, b_dash = interval_update(A, B, c_dash)
        else:
            a_dash, b_dash = A, B
        return a_dash, b_dash

    def bracket(c):
        """Generate initial interval [a, b] that satisfies the opposite
        slope condition (dphi(a) < 0, dphi(b) > 0).
        """
        cs = list()
        for j in range(10):
            cs.append(c)

            dphi_j = get_phi_dphi("g", c)

            if (dphi_j >= 0) and (j == 0):
                return 0, c

            phi_j = get_phi_dphi("f", c)
            if dphi_j >= 0:
                phi_inds = np.array([get_fg("f", c) for c in cs[:-1]]) <= (phi0 + epsk)
                # See https://stackoverflow.com/a/8768734
                ci = len(phi_inds) - phi_inds[::-1].argmax() - 1
                return cs[ci], c
            elif phi_j > (phi0 + epsk):
                return bisect(0, c)

            c *= rho

    def norm_inf(arr):
        """Returns infinity norm of given array."""
        return np.linalg.norm(arr, np.inf)

    def initial():
        """Get an initial guess for alpha."""
        if (~np.isclose(x0, np.zeros_like(x0))).any():
            c = psi_0 * norm_inf(x0)/norm_inf(g0)
        elif not np.isclose(f0, 0):
            c = psi_0 * f0 / norm_inf(g0)**2
        else:
            c = 1
        return c

    def take_quad_step(alpha, g0_):
        """Try to get alpha for minimum step from quadratic interpolation."""
        import pdb; pdb.set_trace()
        fact = max(psi_low, g0_/(dphi0*psi_2))
        alpha_ = min(fact, psi_hi) * alpha
        phi_ = get_phi_dphi("f", alpha_)
        denom = 2*((phi_-phi0)/alpha_ - dphi0)
        f_temp = get_fg("f", alpha_)
        if denom > 0.:
            c = -dphi0*alpha_ / denom
            if f_temp > get_fg("f", 0):
                c = max(c, alpha_*1e-10)
        else:
            c = alpha
        return c

    if alpha_init is None and alpha_prev:
        alpha_init = alpha_prev
    if alpha_init is None and alpha_prev is None:
        alpha_init = initial()

    # Put everything in a try/except block because now everytime
    # we evaluate phi/dphi at some alpha and both phi and dphi
    # are present for this alpha, e.g. from a previous calculation,
    # convergence of the linesearch will be checked, and
    # LineSearchConverged may be raised. Using exceptions enables
    # us to also return from nested functions.
    try:
        if quad_step:
            g0_ = -2*abs(get_fg("f", 0)/alpha_init) if (dphi0_prev is None) \
                  else dphi0_prev
            alpha_init = take_quad_step(psi_2*alpha_init, g0_)
        # This may raise LineSearchConverged
        _ = get_phi_dphi("fg", alpha_init)

        # TODO: cubic interpolation for better alpha_init
        ak, bk = bracket(alpha_init)
        for k in range(max_cycles):
            if cond(ak):
                break
            # secant² step
            a, b = double_secant(ak, bk)
            if (b - a) > gamma*(bk - ak):
                # Bisection step
                c = (a + b)/2
                a, b = interval_update(a, b, c)
            ak, bk = a, b
    except LineSearchConverged as lsc:
        ak = lsc.alpha

    f_new, g_new = get_fg("fg", ak)
    return ak, f_new, g_new, dphi0


@linesearch_wrapper("armijo")
def backtracking(x0, p, get_phi_dphi, get_fg, conds, max_cycles,
                 alpha_init=1., rho_lo=5e-2, rho_hi=0.9):
    """Backtracking line search enforcing Armijo conditions.

    Uses only energy evaluations.

    See [1], Chapter 3, Line Search methods, Section 3.1 p. 31 and
    Section 3.5 p. 56."""
    cond = conds["armijo"]

    log("Starting backtracking line search")
    phi0, dphi0 = get_phi_dphi("fg", 0)

    alpha_prev = None
    alpha = alpha_init
    for i in range(max_cycles):
        phi_i = get_phi_dphi("f", alpha)
        log(f"\tCycle {i:02d}: alpha={alpha:.6f}, ϕ={phi_i:.6f}")

        if cond(alpha):
            log(f"\tLine search converged after {i} cycles.")
            break

        if i == 0:
            # Quadratic interpolation
            alpha_new = interpol_alpha_quad(phi0, dphi0, phi_i, alpha)
            type_ = "Quadratic"
        else:
            # Cubic interpolation
            phi_prev = get_phi_dphi("f", alpha_prev)
            alpha_new = interpol_alpha_cubic(phi0, dphi0, phi_prev, phi_i, alpha_prev, alpha)
            type_ = "Cubic"
        log(f"\tNew alpha from {type_}: {alpha_new:.6f}")

        lower_bound = alpha * rho_lo
        upper_bound = alpha * rho_hi
        if alpha_new < lower_bound:
            log("\tNew alpha is too big!")
        if alpha_new > upper_bound:
            log("\tNew alpha is too high!")

        # Assert that alpha doesn't change too much compared to the previous alpha   
        alpha_new = min(alpha_new, upper_bound)
        alpha_new = max(alpha_new, lower_bound)
        alpha_prev = alpha
        alpha = alpha_new
        log(f"\tAlpha for next cycles: {alpha:.6f}\n")
    else:
        raise LineSearchNotConverged

    return alpha


@linesearch_wrapper("wolfe")
def wolfe(x0, p, get_phi_dphi, get_fg, conds, max_cycles,
          alpha_init=1., alpha_min=0.01, alpha_max=100., fac=2):
    """Wolfe line search.

    Uses only energy & gradient evaluations.

    See [1], Chapter 3, Line Search methods, Section 3.5 p. 60."""

    phi0, dphi0 = get_phi_dphi("fg", 0)

    def zoom(alpha_lo, alpha_hi, phi_lo,
             phi_alpha_=None, alpha_0_=None, max_cycles=10):

        alphas = list()
        phi_alphas = list()
        if phi_alpha_:
            phi_alphas = [phi_alpha_, ]
        if alpha_0_:
            alphas = [alpha_0_, ]

        for j in range(max_cycles):
            # Interpoaltion of alpha between alpha_lo, alpha_hi
            #
            # Try cubic interpolation if at least two additional alphas and
            # corresponding phi_alpha values are available beside alpha = 0.
            if len(phi_alphas) > 1:
                alpha_prev = alphas[-1]
                phi_alpha_prev = phi_alphas[-1]
                alpha_j = interpol_alpha_cubic(phi0, dphi0,
                                               phi_alpha_, phi_alpha_prev,
                                               alpha_0_, alpha_prev
                )
            # Try quadratic interpolation if at one additional alpha and
            # corresponding phi_alpha value is available beside alpha = 0.
            elif len(phi_alphas) == 1:
                alpha_j = interpol_alpha_quad(phi0, dphi0, phi_alpha_, alpha_0_)
            # Fallback to simple bisection
            else:
                alpha_j = (alpha_lo + alpha_hi) / 2

            phi_j = get_phi_dphi("f", alpha_j)
            # Store the values so they can be reused for cubic interpolation
            alphas.append(alpha_j)
            phi_alphas.append(phi_j)

            # True if alpha is still too big or if the function value
            # increased compared to the previous cycle.
            if (not conds["armijo"](alpha_j) or phi_j > phi_lo):
                # Shrink interval to (alpha_lo, alpha_j)
                alpha_hi = alpha_j
                continue

            dphi_j = get_phi_dphi("g", alpha_j)
            if conds["curv"](alpha_j):
                print(f"\tzoom converged after {j+1} cycles.")
                return alpha_j

            if (dphi_j * (alpha_hi - alpha_lo)) >= 0:
                alpha_hi = alpha_lo
            # Shrink interval to (alpha_j, alpha_hi)
            alpha_lo = alpha_j
        raise Exception("zoom() didn't converge in {j+1} cycles!")

    alpha_prev = 0
    phi_prev = phi0
    if alpha_init is not None:
        alpha_i = alpha_init
    # This does not seem to help
    # elif f_0_prev is not None:
        # alpha_i = min(1.01*2*(f_0 - f_0_prev) / dphi_0, 1.)
        # print("ai", alpha_i)
        # alpha_i = 1. if alpha_i < 0. else alpha_i
    else:
        alpha_i = 1.0

    try:
        for i in range(10):
            phi_i = get_phi_dphi("f", alpha_i)
            if (not conds["armijo"](alpha_i) or ((phi_i >= phi_prev) and i > 0)):
                zoom(alpha_prev, alpha_i, phi_prev, phi_i, alpha_i)

            dphi_i = get_phi_dphi("g", alpha_i)
            if conds["curve"](alpha_i):
                raise LineSearchConverged(alpha_i)

            if dphi_i >= 0:
                zoom(alpha_i, alpha_prev, phi_i, phi_alpha_=phi_i, alpha_0_=alpha_i)
            prev_alpha = alpha_i
            alpha_i = min(fac * alpha_i, alpha_max)
        else:
            raise LineSearchNotConverged
    except LineSearchConverged as lsc:
        alpha = lsc.alpha

    return alpha
