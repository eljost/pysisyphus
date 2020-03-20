import matplotlib.pyplot as plt
import numpy as np
import pytest


from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Rosenbrock import Rosenbrock
from pysisyphus.optimizers.line_searches import hager_zhang, backtracking


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
            # "g0": g,
            # "alpha_init": 0.25,
        # }
        # alpha, f_new, g_new = backtracking_linesearch(**bt_kwargs)
        # # alpha, f_new, g_new = btd_(**bt_kwargs)

        kwargs = {
            "f": fun,
            "df": jac,
            "x0": x,
            "p": d,
            "f0": f,
            # "df0": g,
            # "alpha_init": None,
            # "alpha_init": 0.230871 if alpha_prev is None else None,
            "alpha_prev": alpha_prev,
            "f_prev": f_prev,
            "quad_step": True,
            # dphi0_prev will be set if alpha_prev is not None
            "dphi0_prev": None if alpha_prev is None else dphi0_prev,  # noqa: F821
        }
        # alpha, f_new, g_new, dphi0_prev = hager_zhang_linesearch(**kwargs)
        kwargs.update({"g0": g})
        alpha, f_new, g_new, dphi0_prev = hager_zhang(**kwargs)

        x_new = x + alpha*d

        y = g_new - g
        dy = d.dot(y)
        beta = y.dot(g_new)/dy - y.dot(y)/dy * d.dot(g_new)/dy

        eta = 0.4
        etak = eta * d.dot(g)/d.dot(d)
        beta = max(beta, etak)
        # print(f"\tbeta={beta:.6f}")
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
def test_cg_descent(calc, x0):
    x0 = np.array(x0)
    geom = calc.get_geom(x0)

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

    xs = cg_descent(x0, fun, jac)#, max_cycles=1)
    print(f"func calls: {fun.calls}")
    print(f"jac calls: {jac.calls}")
    g_inf = np.linalg.norm(jac(xs[-1]), np.inf)
    print(f"f_min={fun(xs[-1]):.8e}, ||g_min||_inf={g_inf:.8e}, x_min={xs[-1]}")

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
    test_cg_descent()
