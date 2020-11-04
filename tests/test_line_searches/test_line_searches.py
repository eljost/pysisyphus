import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Rosenbrock import Rosenbrock
from pysisyphus.calculators.AnaPotBase import AnaPotBase
from pysisyphus.line_searches import Backtracking, HagerZhang, StrongWolfe
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.PreconSteepestDescent import PreconSteepestDescent
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


class CGDescent(Optimizer):
    def __init__(self, geometry, alpha_init=0.5, **kwargs):
        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init

        self.alpha_prev = None

    def prepare_opt(self):
        forces = self.geometry.forces
        self.step_direction = forces / np.linalg.norm(forces)

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        f = lambda coords: self.geometry.get_energy_at(coords)
        df = lambda coords: -self.geometry.get_energy_and_forces_at(coords)["forces"]

        try:
            f_prev = self.energies[-2]
        except IndexError:
            f_prev = None

        alpha_init = self.alpha_init if self.alpha_prev is None else None

        # kwargs = {
        # "f": f,
        # "df": df,
        # "x0": self.geometry.coords,
        # "p": self.step_direction,
        # "f0": energy,
        # "g0": -forces,
        # "alpha_init": alpha_init,
        # "alpha_prev": self.alpha_prev,
        # "f_prev": f_prev,
        # "quad_step": True,
        # # dphi0_prev will be set if alpha_prev is not None
        # "dphi0_prev": None if not self.alpha_prev else self.dphi0_prev,  # noqa: F821
        # }
        # alpha, f_new, g_new, dphi0_prev = hager_zhang(**kwargs)

        kwargs = {
            "geometry": self.geometry,
            "p": self.step_direction,
            "f0": energy,
            "g0": -forces,
            "alpha_init": alpha_init,
            "alpha_prev": self.alpha_prev,
            "f_prev": f_prev,
            "quad_step": True,
            # dphi0_prev will be set if alpha_prev is not None
            "dphi0_prev": None
            if not self.alpha_prev
            else self.dphi0_prev,  # noqa: F821
        }
        line_search = HagerZhang(**kwargs)
        hz_result = line_search.run()
        alpha = hz_result.alpha
        f_new = hz_result.f_new
        g_new = hz_result.g_new
        dphi0_prev = hz_result.dphi0

        self.dphi0_prev = dphi0_prev
        step = alpha * self.step_direction

        d = self.step_direction
        y = g_new - -forces
        dy = d.dot(y)
        beta = y.dot(g_new) / dy - y.dot(y) / dy * d.dot(g_new) / dy

        eta = 0.4
        etak = eta * d.dot(-forces) / d.dot(d)
        beta = max(beta, etak)
        # print(f"\tbeta={beta:.6f}")
        self.step_direction = -g_new + beta * d

        self.alpha_prev = alpha

        return step


@pytest.mark.parametrize(
    "line_search, ref_cycle",
    [
        # ("armijo", 40),
        ("strong_wolfe", 13),
        ("hz", 9),
        (None, 18),
    ],
)
def test_line_search(line_search, ref_cycle):
    geom = AnaPot.get_geom((0.687, 1.57, 0.0))

    opt_kwargs = {
        "thresh": "gau_tight",
        "line_search": line_search,
        "max_cycles": 64,
        "precon": False,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    # geom.calculator.plot_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(0.98555442)


def test_1d_steepest_descent():
    V_str = "2*x**4 + 5*x**3 - 2*x**2 + 10*x"
    geom = AnaPotBase.get_geom((-3, 0.0, 0.0), V_str=V_str)

    opt = PreconSteepestDescent(geom, thresh="gau_tight", precon=False)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 16


@pytest.mark.parametrize(
    "calc, x0, ref_cycle",
    [
        (AnaPot, (0, 3.0, 0), 10),
        (Rosenbrock, (-1.2, 1.0, 0.0), 34),
    ],
)
def test_cg_descent(calc, x0, ref_cycle):
    x0 = np.array(x0)
    geom = calc.get_geom(x0)
    calc = geom.calculator

    opt_kwargs = {
        "thresh": "gau_tight",
    }
    opt = CGDescent(geom, **opt_kwargs)
    opt.run()

    # geom.calculator.plot_opt(opt)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle


@pytest.mark.parametrize(
    # "use_grad, ref_cycle", [
    "line_search, ref_cycle",
    [
        ("armijo", 7),
        ("armijo_fg", 6),
    ],
)
def test_armijo(line_search, ref_cycle):
    geom = AnaPot.get_geom((0.2, 1.3, 0.0))

    opt_kwargs = {
        "thresh": "gau_tight",
        "precon": False,
        # "use_grad": use_grad,
        "line_search": line_search,
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    # geom.calculator.plot_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(-0.51340926)
