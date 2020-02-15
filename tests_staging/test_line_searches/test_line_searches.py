import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.Rosenbrock import Rosenbrock
from pysisyphus.calculators.AnaPotBase import AnaPotBase
from pysisyphus.optimizers.line_searches import backtracking, wolfe, hager_zhang
from pysisyphus.optimizers.Optimizer import Optimizer

from pysisyphus.line_searches import Backtracking, HagerZhang, StrongWolfe


class SteepestDescent(Optimizer):

    def __init__(self, geometry, alpha_init=0.5, line_search="armijo",
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init
        self.line_search = line_search

        ls_cls = {
            "armijo": Backtracking,
            "strong_wolfe": StrongWolfe,
            "hz": HagerZhang,
        }
        # ls_funcs = {
            # "armijo": backtracking,
            # "strong_wolfe": wolfe,
        # }
        # self.line_search_func = ls_funcs[self.line_search]
        self.line_search_cls = ls_cls[self.line_search]

        self.alpha_prev = None

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        step_dir = forces / np.linalg.norm(forces)

        # f = lambda coords: self.geometry.get_energy_at(coords)
        # df = lambda coords: self.geometry.get_energy_and_forces_at(coords)["forces"]
        # kwargs = {
            # "f": f,
            # "df": df,
            # "x0": self.geometry.coords,
            # "p": step_dir,
            # "f0": energy,
            # "g0": -forces,
            # # "alpha_init": self.alpha_prev if self.alpha_prev else self.alpha_init,
            # "alpha_init": self.alpha_init,
        # }
        # alpha = self.line_search_func(**kwargs)

        # OO Interface
        kwargs = {
            "geometry": self.geometry,
            "p": step_dir,
            "f0": energy,
            "g0": -forces,
            "alpha_init": self.alpha_init,
        }
        line_search = self.line_search_cls(**kwargs)
        line_search_result = line_search.run()
        alpha = line_search_result.alpha

        step = alpha * step_dir
        self.alpha_prev = alpha
        return step


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
            "dphi0_prev": None if not self.alpha_prev else self.dphi0_prev,  # noqa: F821
        }
        line_search = HagerZhang(**kwargs)
        hz_result = line_search.run()
        alpha = hz_result.alpha
        f_new = hz_result.f_new
        g_new = hz_result.g_new
        dphi0_prev = hz_result.dphi0

        self.dphi0_prev = dphi0_prev
        step = alpha*self.step_direction

        d = self.step_direction
        y = g_new - -forces
        dy = d.dot(y)
        beta = y.dot(g_new)/dy - y.dot(y)/dy * d.dot(g_new)/dy

        eta = 0.4
        etak = eta * d.dot(-forces)/d.dot(d)
        beta = max(beta, etak)
        # print(f"\tbeta={beta:.6f}")
        self.step_direction = -g_new + beta*d

        self.alpha_prev = alpha

        return step


@pytest.mark.parametrize(
    "line_search, ref_cycle, ref_energy",
    [
        ("armijo", 32, 0.98555442),
        ("strong_wolfe", 57, 0.98555442),
        ("hz", 63, 0.98555442),
    ]
)
def test_line_search(line_search, ref_cycle, ref_energy):
    geom = AnaPot.get_geom((0.687, 1.57, 0.))

    opt_kwargs = {
        "thresh": "gau_tight",
        "line_search": line_search,
        "max_cycles": 75,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()

    # geom.calculator.plot_opt(opt)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(ref_energy)


def test_1d_steepest_descent():
    V_str = "2*x**4 + 5*x**3 - 2*x**2 + 10*x"
    geom = AnaPotBase.get_geom((-3, 0., 0.), V_str=V_str)

    opt = SteepestDescent(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 5


@pytest.mark.parametrize(
    "calc, x0, ref_cycle",
    [
        (AnaPot, (0, 3., 0), 10),
        (Rosenbrock, (-1.2, 1.0, 0.), 34),
])
def test_cg_descent(calc, x0, ref_cycle):
    x0 = np.array(x0)
    geom = calc.get_geom(x0)
    calc = geom.calculator

    opt_kwargs =  {
        "thresh": "gau_tight",
    }
    opt = CGDescent(geom, **opt_kwargs)
    opt.run()

    # geom.calculator.plot_opt(opt)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
