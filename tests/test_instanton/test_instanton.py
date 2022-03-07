import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.helpers import geom_loader
from pysisyphus.irc.Instanton import Instanton
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer


def fin_diff(inst, fun_name, key, deriv_arr, dr=1e-6):
    dr = 1e-6
    coords = inst.coords.copy()
    num_grad = np.zeros_like(coords)
    fun = getattr(inst, fun_name)

    def calc():
        return fun()[key]

    for i, crds in enumerate(coords):
        minus = coords.copy()
        minus[i] -= dr
        inst.coords = minus
        S_minus = calc()
        plus = coords.copy()
        plus[i] += dr
        inst.coords = plus
        S_plus = calc()
        dS = (S_plus - S_minus) / (2 * dr)
        deriv_arr[i] = dS
    # Restore original coordinates
    inst.coords = coords.copy()
    return deriv_arr


def test_instanton_action():
    geom = AnaPot().get_saddles(i=0, geom_kwargs={"coord_type": "mwcartesian"})
    calc = geom.calculator

    def calc_getter():
        return AnaPot()

    P = 10
    inst = Instanton.from_ts(geom, calc_getter=calc_getter, P=P)

    res = inst.action_gradient()
    action = res["action"]
    assert action == pytest.approx(12.09601)

    grad_res = inst.action_gradient()
    gradient = grad_res["gradient"]

    dr = 1e-6
    coords = inst.coords.copy()
    # Check gradient
    num_grad = np.zeros_like(coords)
    fin_diff(inst, "action", "action", num_grad)
    np.testing.assert_allclose(gradient, num_grad)

    hess_res = inst.action_hessian()
    hessian = hess_res["hessian"]
    # Check Hessian
    num_hessian = np.zeros((coords.size, coords.size))
    fin_diff(inst, "action_gradient", "gradient", num_hessian)
    np.testing.assert_allclose(hessian, num_hessian, atol=1e-8)


def test_instanton_opt():
    scale = 1e-2
    geom = AnaPot().get_saddles(
        i=0,
        geom_kwargs={
            "coord_type": "mwcartesian",
        },
        calc_kwargs={
            "scale": scale,
        },
    )
    calc = geom.calculator

    def calc_getter():
        return AnaPot(scale=scale)

    P = 20
    T = 100
    inst = Instanton.from_ts(geom, P=P, calc_getter=calc_getter, T=T)

    opt = RSIRFOptimizer(inst, hessian_init="calc", hessian_recalc=5)
    opt.run()

    # coords = np.array(opt.coords)
    # coords = coords.reshape(-1, P, 3)
    # calc.plot()
    # ax = calc.ax
    # ax.plot(*coords[0].T[:2], "o-", label="start")
    # ax.plot(*coords[-1].T[:2], "o-", label="end")
    # ax.legend()
    # plt.show()

    assert opt.cur_cycle == 21
    assert inst.energy == pytest.approx(1.66692135)


def test_sequential_cooling():
    calc_kwargs = {"scale": 1e-2}
    geom = AnaPot().get_saddles(
        i=0,
        geom_kwargs={
            "coord_type": "mwcartesian",
        },
        calc_kwargs=calc_kwargs,
    )
    calc = geom.calculator

    def calc_getter():
        return AnaPot(**calc_kwargs)

    P = 20
    T = 400
    Ts = np.linspace(T, 100, num=10)
    inst = Instanton.from_ts(geom, P=P, calc_getter=calc_getter)
    init = list()
    fin = list()
    for i, T in enumerate(Ts):
        print(f"@@@ CYCLE {i}, T={T} K")
        inst = Instanton.from_instanton(inst, calc_getter=calc_getter, T=T)
        init.append(inst.coords.copy())
        opt = RSIRFOptimizer(
            inst, hessian_init="calc", hessian_recalc=5
        )
        opt.run()
        fin.append(inst.coords.copy())

    # init = np.array(inst)
    # fin = np.array(fin)
    # calc.plot()
    # ax = calc.ax
    # for i, (coords, T) in enumerate(zip(fin, Ts)):
    # c3d = coords.reshape(-1, 3)
    # ax.plot(*c3d.T[:2], "o-", label=f"{T:.1f} K")
    # ax.legend()
    # plt.show()
