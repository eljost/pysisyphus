import matplotlib.pyplot as plt
import numpy as np
import pytest
from tensorflow_probability.python.optimizer.linesearch import hager_zhang as tfhz


from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Rosenbrock import Rosenbrock
from pysisyphus.optimizers.line_search import interpol_alpha_cubic, interpol_alpha_quad


class LinesearchConverged(Exception):

    def __init__(self, alpha):
        self.alpha = alpha


def backtracking_linesearch(f, df, x0, p, f0=None, df0=None, alpha_init=None,
                            c1=1e-4, rho_lo=0.1, rho_hi=0.5,):
    alpha_fs = {}
    def phi(alpha):
        try:
            f_alpha = alpha_fs[alpha]
        except KeyError:
            f_alpha = f(x0 + alpha*p)
            alpha_fs[alpha] = f_alpha
        return f_alpha

    alpha_dfs = {}
    dphis = {}
    def dphi(alpha):
        try:
            df_alpha = alpha_dfs[alpha]
            dphi_ = df_alpha @ p
        except KeyError:
            df_alpha = df(x0 + alpha*p)
            alpha_dfs[alpha] = df_alpha
            dphi_ = df_alpha @ p
            dphis[alpha] = dphi_
        return dphi_

    if f0 is None:
        phi0 = phi(0)
    else:
        phi0 = f0
        alpha_fs[0.] = f0
    if df0 is None:
        dphi0 = dphi(0)
    else:
        dphi0 = df0 @ p
        alpha_dfs[0.] = df0

    def sufficiently_decreased(alpha):
        return phi(alpha) <= (phi0 + c1 * alpha * dphi0)

    alpha = alpha_init
    alpha_prev = alpha_init
    for i in range(25):
        phi_i = phi(alpha)

        if sufficiently_decreased(alpha):
            print(f"\tbacktracking converged after {i} cycles.")
            break

        if i == 0:
            # Quadratic interpolation
            alpha_new = interpol_alpha_quad(phi0, dphi0, phi_i, alpha)
            print(f"\tQUAD alpha={alpha:.6f}")
        else:
            # Cubic interpolation
            phi_prev = phi(alpha_prev)
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


    dphi_i = dphi(alpha)
    return alpha, alpha_fs[alpha], alpha_dfs[alpha]


def hager_zhang_linesearch(f, df, x0, p, f0=None, df0=None, alpha_init=None,
                           alpha_prev=None, f_prev=None, dphi0_prev=None, quad_step=False,
                           c1=.1, c2=0.9, eps=1e-6, theta=0.5, gamma=0.5, rho=5,
                           psi_0=.01, psi_1=.1, psi_2=2., psi_low=0.1, psi_hi=10,
                           Delta=.7, omega=1e-3, max_cycles=10, max_bisects=10):
    alpha_fs = {}
    def phi(alpha):
        try:
            f_alpha = alpha_fs[alpha]
        except KeyError:
            f_alpha = f(x0 + alpha*p)
            alpha_fs[alpha] = f_alpha
        return f_alpha

    alpha_dfs = {}
    dphis = {}
    def dphi(alpha):
        try:
            df_alpha = alpha_dfs[alpha]
            dphi_ = df_alpha @ p
        except KeyError:
            df_alpha = df(x0 + alpha*p)
            alpha_dfs[alpha] = df_alpha
            dphi_ = df_alpha @ p
            dphis[alpha] = dphi_
        return dphi_


    def get_phi_dphi(what, alpha):
        whats = "f g fg".split()
        assert what in whats
        calc_funcs = {
            "f": phi,
            "g": dphi,
        }
        result = [calc_funcs[w](alpha) for w in what]
        # Check if we got both phi and dphi for alpha now. If so we
        # can check if the chosen condition (Wolfe/approx. Wolfe) is
        # satisfied.
        if (alpha in alpha_fs) and (alpha in alpha_dfs) and cond(alpha):
            raise LinesearchConverged(alpha)
        # Dont return a list if only f or g was requested.
        if len(what) == 1:
            result = result[0]
        return result

    if f0 is None:
        phi0 = phi(0)
    else:
        phi0 = f0
        alpha_fs[0.] = f0
    if df0 is None:
        dphi0 = dphi(0)
    else:
        dphi0 = df0 @ p
        alpha_dfs[0.] = df0

    epsk = eps * abs(alpha_fs[0.])
    # wolfe_hi = c1*dphi0
    # wolfe_low = c2*dphi0

    def sufficiently_decreased(alpha):
        return phi(alpha) <= (phi0 + c1 * alpha * dphi0)

    def curvature_condition(alpha):
        return dphi(alpha) >= c2*dphi0

    def wolfe_condition(alpha):
        return sufficiently_decreased(alpha) and curvature_condition(alpha)

    def approx_wolfe_condition(alpha):
        return (2*c1-1)*dphi0 >= dphi(alpha) >= c2*dphi0

    def t2_condition(alpha):
        phi(alpha) <= phi0 + epsk

    def bisect(a, b):
        for i in range(max_bisects):
            # U3 a.
            d = (1 - theta)*a + theta*b
            # dphi_d = dphi(d)
            dphi_d = get_phi_dphi("g", d)
            if dphi_d >= 0:
                return a, d

            # phi_d = phi(d)
            phi_d = get_phi_dphi("f", d)
            # U3 b.
            if (dphi_d < 0) and (phi_d <= phi0 + epsk):
                a = d
            # U3 c.
            elif (dphi_d < 0) and (phi_d > phi0 + epsk):
                b = d
        raise Exception("Bisect failed!")

    def interval_update(a, b, c):
        # U0
        if not (a < c < b):
            return a, b

        # phi_c = phi(c)
        # dphi_c = dphi(c)
        phi_c, dphi_c = get_phi_dphi("fg", c)
        # U1, sign of slope projection changed. We already passed the minimum.
        if dphi_c >= 0:
            return a, c
        # U2, we are moving towards the minimum.
        elif (dphi_c < 0) and (phi_c <= phi0 + epsk):
            return c, b

        # U3, phi_c increased above phi0, so we probably passed the minimum.
        return bisect(a, c)

    def secant(a, b):
        return (a*dphi(b) - b*dphi(a)) / (dphi(b) - dphi(a))

    def double_secant(a, b):
        fc = f.calls
        dfc = df.calls
        c = secant(a, b)
        fc_ = f.calls
        dfc_ = df.calls
        A, B = interval_update(a, b, c)
        fc__ = f.calls
        dfc__ = df.calls
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

    # def pri(a, b):
        # pa = phi(a)
        # pb = phi(b)
        # fa = alpha_fs[a]
        # fb = alpha_fs[b]
        # dpa = dphi(a)
        # dpb = dphi(b)
        # da = dpa-wolfe_hi
        # db = dpb-wolfe_hi
        # print(f"a={a:>10.6f} b={b:>10.6f} fa={fa:>10.6f} fb={fb:>10.6f} "
              # f"da={da:>10.6f} db={db:>10.6f}"
        # )

    def bracket(c):
        cs = list()
        for j in range(10):
            cs.append(c)

            dphi_j = get_phi_dphi("g", c)
            # dphi_j = dphi(c)

            if (dphi_j >= 0) and (j == 0):
                # pri(0, c)
                return 0, c

            phi_j = get_phi_dphi("f", c)
            # phi_j = phi(c)
            if dphi_j >= 0:
                phi_inds = np.array([alpha_fs[c] for c in cs[:-1]]) <= (phi0 + epsk)
                # See https://stackoverflow.com/a/8768734
                ci = len(phi_inds) - phi_inds[::-1].argmax() - 1
                # pri(cs[ci], c)
                return cs[ci], c
            elif (dphi_j < 0) and (phi_j > (phi0 + epsk)):
                # pri(*bisect(0, c))
                return bisect(0, c)

            c *= rho

    def norm_inf(arr): return np.linalg.norm(arr, np.inf)

    def initial():
        if (~np.isclose(x0, np.zeros_like(x0))).any():
            c = psi_0 * norm_inf(x0)/norm_inf(df0)
        elif not np.isclose(f0, 0):
            c = psi_0 * f0 / norm_inf(df0)**2
        else:
            c = 1
        print(f"I0 alpha={c:.6f}")
        return c

    def take_quad_step(alpha, df0_):
        fact = max(psi_low, df0_/(dphi0*psi_2))
        print(f"t={fact:.6f}")
        alpha_ = min(fact, psi_hi) * alpha
        print(f"alpha for quadstep={alpha_:.6f}")
        # alpha_ = psi_1 * 2 * alpha
        phi_ = get_phi_dphi("f", alpha_)
        num = dphi0*alpha_**2
        # denom = 2*(phi_ - phi0 - dphi0*alpha_)
        denom = 2*((phi_-phi0)/alpha_ - dphi0)
        f_temp = alpha_fs[alpha_]
        print(f"f_temp={f_temp:.6f}")
        if denom > 0.:
            c_ = -dphi0*alpha_ / denom
            if f_temp > alpha_fs[0]:
                c_ = max(c_, alpha_*1e-10)
        else:
            c_ = alpha
        return c_

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
    cond = wolfe_condition
    # cond = t2_condition
    # cond = approx_wolfe_condition

    if alpha_init is None and alpha_prev:
        alpha_init = alpha_prev
    if alpha_init is None and alpha_prev is None:
        alpha_init = initial()

    # if cond(alpha_init):
        # ak = alpha_init
        # print(f"\tfinal quad alpha={ak:.6f}")
        # return ak, alpha_fs[ak], alpha_dfs[ak]

    if quad_step:
        df0 = -2*abs(alpha_fs[0]/alpha_init) if (dphi0_prev is None) else dphi0_prev
        alpha_init = take_quad_step(psi_2*alpha_init, df0)
        # print("start ", end=""); pri(0, alpha_init)

        if cond(alpha_init):
            ak = alpha_init
            print(f"\tfinal quad alpha={ak:.6f}")
            return ak, alpha_fs[ak], alpha_dfs[ak], dphi0


    print(f"\talpha_init={alpha_init:.6f}")
    # Put everything in a try/except block because now everytime
    # we evaluate phi/dphi at some alpha and both phi and dphi
    # are present convergence of the linesearch will be checked,
    # and LinesearchConverged may be raised.
    try:
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
    
    print(f"\tfinal alpha={ak:.6f}")
    return ak, alpha_fs[ak], alpha_dfs[ak], dphi0


def sd(x0, fun, jac, max_cycles=50):
    f0 = fun(x0)
    g0 = jac(x0)
    p0 = -g0
    
    x = x0
    f = f0
    g = g0
    p = p0
    
    fs = list()
    xs = list()
    for i in range(max_cycles):
        xs.append(x)
        fs.append(f)
        
        inf_norm = np.linalg.norm(g, np.inf)
        print(f"Cycle {i:02d}: {inf_norm:.6f}")
        if inf_norm <= 1e-8:
            print("Converged")
            break

        bt_kwargs = {
            "f": fun,
            "df": jac,
            "x0": x,
            "p": p,
            "f0": f,
            "df0": g,
            "alpha_init": 0.25,
        }
        alpha, f_new, g_new = backtracking_linesearch(**bt_kwargs)
        x_new = x + alpha*p
        p_new = -g_new

        x = x_new
        f = f_new
        g = g_new
        p = p_new
    return np.array(xs)


def cg_descent(x0, fun, jac, max_cycles=50):
    x = x0
    f = fun(x0)
    g = jac(x0)
    d = -g
    
    xs = list()
    alpha_prev = None
    f_prev = None
    for i in range(max_cycles):
        xs.append(x)
        
        inf_norm = np.linalg.norm(g, np.inf)
        print(f"Cycle {i:02d}: {inf_norm:.6f}")
        if inf_norm <= 1e-8:
            print("Converged")
            break

        # bt_kwargs = {
            # "f": fun,
            # "df": jac,
            # "x0": x,
            # "p": d,
            # "f0": f,
            # "df0": g,
            # "alpha_init": 0.25,
        # }
        # alpha, f_new, g_new = backtracking_linesearch(**bt_kwargs)

        kwargs = {
            "f": fun,
            "df": jac,
            "x0": x,
            "p": d,
            "f0": f,
            "df0": g,
            # "alpha_init": None,
            # "alpha_init": 0.230871 if alpha_prev is None else None,
            "alpha_prev": alpha_prev,
            "f_prev": f_prev,
            "quad_step": True,
            "dphi0_prev": None if alpha_prev is None else dphi0_prev,
        }
        alpha, f_new, g_new, dphi0_prev = hager_zhang_linesearch(**kwargs)

        x_new = x + alpha*d

        y = g_new - g
        dy = d.dot(y)
        beta = y.dot(g_new)/dy - y.dot(y)/dy * d.dot(g_new)/dy

        eta = 0.4
        etak = eta * d.dot(g)/d.dot(d)
        beta = max(beta, etak)
        print(f"\tbeta={beta:.6f}")
        d_new = -g_new + beta*d

        f_prev = f
        alpha_prev = alpha

        x = x_new
        f = f_new
        g = g_new
        d = d_new
    else:
        # Append last coords if max_cycles was reached
        xs.append(x)
    return np.array(xs)


@pytest.mark.parametrize(
    "calc, x0",
    [
        pytest.param(AnaPot, (0, 3., 0)),
        pytest.param(Rosenbrock, (-1.2, 1.0, 0.)),
])
def test_hager_zhang(calc, x0):
    # x0 = np.array((0, 3, 0.))
    x0 = np.array(x0)
    # geom = AnaPot.get_geom(x0)
    geom = calc.get_geom(x0)
    # x0 = np.array((-1.2, 1.0, 0.))
    # geom = Rosenbrock.get_geom(x0)

    calc = geom.calculator

    def fun(x0):
        if not hasattr(fun, "calls"):
            fun.calls = 0
        res = geom.get_energy_and_forces_at(x0)
        fun.calls += 1
        return res["energy"]


    def jac(x0):
        if not hasattr(jac, "calls"):
            jac.calls = 0
        res = geom.get_energy_and_forces_at(x0)
        jac.calls += 1
        return -res["forces"]

    # xs = sd(x0, fun, jac)
    xs = cg_descent(x0, fun, jac)#, max_cycles=2)
    print(f"func calls: {fun.calls}")
    print(f"jac calls: {jac.calls}")
    g_inf = np.linalg.norm(jac(xs[-1]), np.inf)
    print(f"f_min={fun(xs[-1]):.8e}, ||g_min||_inf={g_inf:.8e}, x_min={xs[-1]}")
    print("xs", xs)

    # calc.plot()
    # ax = calc.ax
    # ax.plot(*xs.T[:2], "o-")
    # ax.set_xlim(-1.25, 1)
    # ax.set_ylim(0.5, 3)
    # plt.show()

    # calc.plot()
    # ax = calc.ax
    # ax.plot(*xs.T[:2], "o-")
    # ax.set_xlim(-2.5, 2.5)
    # ax.set_ylim(-0.5, 1.5)
    # plt.show()


if __name__ == "__main__":
    test_hager_zhang()
