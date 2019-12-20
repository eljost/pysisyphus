import numpy as np

from pysisyphus.optimizers.line_search import interpol_alpha_cubic, interpol_alpha_quad


class LinesearchConverged(Exception):

    def __init__(self, alpha):
        self.alpha = alpha


class LinesearchNotConverged(Exception):
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
                   and (alpha in alpha_gs) and (alpha in alpha_gs) and cond_func(alpha):
                    raise LinesearchConverged(alpha)
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
                "wolfe": wolfe_condition,
                "strong_wolfe": strong_wolfe_condition,
            }

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

            cond_func = conds[cond]

            linesearch_result = func(x0, p, get_phi_dphi, get_fg, cond_func,
                                     max_cycles, *args, **kwargs)
            return linesearch_result
        return wrapper
    return cond_wrapper


@linesearch_wrapper(cond="wolfe")
def hager_zhang(x0, p, get_phi_dphi, get_fg, cond, max_cycles,
                alpha_init=None, alpha_prev=None,
                f_prev=None, dphi0_prev=None, quad_step=False,
                eps=1e-6, theta=0.5, gamma=0.5, rho=5,
                psi_0=.01, psi_1=.1, psi_2=2., psi_low=0.1, psi_hi=10,
                Delta=.7, omega=1e-3, max_bisects=10):
    epsk = eps * abs(get_fg("f", 0.))
    phi0, dphi0 = get_phi_dphi("fg", 0.)
    f0, g0 = get_fg("fg", 0.)

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
    # LinesearchConverged may be raised. Using exceptions enables
    # us to also return from nested functions.
    try:
        if quad_step:
            g0_ = -2*abs(get_fg("f", 0)/alpha_init) if (dphi0_prev is None) \
                  else dphi0_prev
            alpha_init = take_quad_step(psi_2*alpha_init, g0_)
        # This may raise LinesearchConverged
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
    except LinesearchConverged as lsc:
        ak = lsc.alpha

    f_new, g_new = get_fg("fg", ak)
    return ak, f_new, g_new, dphi0


@linesearch_wrapper(cond="armijo")
def backtracking(x0, p, get_phi_dphi, get_fg, cond, max_cycles,
                 alpha_init=None, rho_lo=0.1, rho_hi=0.5):
    phi0, dphi0 = get_phi_dphi("fg", 0)

    alpha = alpha_init
    alpha_prev = alpha_init
    for i in range(max_cycles):
        phi_i = get_phi_dphi("f", alpha)

        if cond(alpha):
            print(f"\tbacktracking converged after {i} cycles.")
            break

        if i == 0:
            # Quadratic interpolation
            alpha_new = interpol_alpha_quad(phi0, dphi0, phi_i, alpha)
            print(f"\tQUAD alpha={alpha:.6f}")
        else:
            # Cubic interpolation
            phi_prev = get_phi_dphi("f", alpha_prev)
            alpha_new = interpol_alpha_cubic(phi0, dphi0, phi_prev, phi_i, alpha_prev, alpha)
            print(f"\tCUB alpha={alpha:.6f}")

        alpha_prev = alpha

        if alpha_new > alpha_prev*rho_hi:
            print("\tProposed alpha is too high!")
        if alpha_new < alpha_prev*rho_lo:
            print("\tProposed alpha is too small!")

        # Assert that alpha doesn't change too much compared to the previous alpha   
        alpha_new = min(alpha_new, alpha_prev*rho_hi)
        alpha = max(alpha_new, alpha_prev*rho_lo)
    else:
        raise LinesearchNotConverged

    # Call this so get_fg will always return something...
    dphi_i = get_phi_dphi("g", alpha, check=False)  # lgtm [py/unused-local-variable]
    f_new, g_new = get_fg("fg", alpha)
    return alpha, f_new, g_new
