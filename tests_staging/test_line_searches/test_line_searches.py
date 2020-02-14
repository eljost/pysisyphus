import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.AnaPotBase import AnaPotBase
from pysisyphus.optimizers.line_searches import backtracking, wolfe, hager_zhang
from pysisyphus.optimizers.Optimizer import Optimizer


class SteepestDescent(Optimizer):

    def __init__(self, geometry, alpha_init=0.5, line_search="armijo",
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init
        self.line_search = line_search

        ls_funcs = {
            "armijo": backtracking,
            "wolfe": wolfe,
            "hz": hager_zhang,
        }
        self.line_search_func = ls_funcs[self.line_search]

        self.alpha_prev = None

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        step_dir = forces / np.linalg.norm(forces)

        f = lambda coords: self.geometry.get_energy_at(coords)
        df = lambda coords: self.geometry.get_energy_and_forces_at(coords)["forces"]

        kwargs = {
            "f": f,
            "df": df,
            "x0": self.geometry.coords,
            "p": step_dir,
            "f0": energy,
            "g0": -forces,
            # "alpha_init": self.alpha_prev if self.alpha_prev else self.alpha_init,
            "alpha_init": self.alpha_init,
        }
        alpha = self.line_search_func(**kwargs)

        step = alpha * step_dir
        self.alpha_prev = alpha
        return step


@pytest.mark.parametrize(
    "line_search, ref_cycle",
    [
        ("armijo", 32),
        # ("wolfe", 32),
        # ("hz", 32),
    ]
)
def test_line_search(line_search, ref_cycle):
    geom = AnaPot.get_geom((0.687, 1.57, 0.))

    opt_kwargs = {
        "thresh": "gau_tight",
        "line_search": line_search,
    }
    opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()

    cs = np.array(opt.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*cs.T[:2])
    plt.show()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

def test_1d_steepest_descent_wolfe():
    V_str = "2*x**4 + 5*x**3 - 2*x**2 + 10*x"
    geom = AnaPotBase.get_geom((-3, 0., 0.), V_str=V_str)

    opt = SteepestDescent(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 5
