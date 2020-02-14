import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.line_searches import backtracking
from pysisyphus.optimizers.Optimizer import Optimizer


class SteepestDescent(Optimizer):

    def __init__(self, geometry, alpha=0.1, **kwargs):
        super().__init__(geometry, **kwargs)

        self.alpha = alpha

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        step_dir = forces / np.linalg.norm(forces)
        # alpha, _, _ = 0.25, None, None

        f = lambda coords: self.geometry.get_energy_at(coords)
        df = lambda coords: None

        kwargs = {
            "f": f,
            "df": df,
            "x0": self.geometry.coords,
            "p": step_dir,
            "f0": energy,
            "g0": -forces,
            "alpha_init": 0.5,
        }
        # alpha, f_new, g_new = backtracking(**kwargs)
        alpha = backtracking(**kwargs)

        step = alpha * step_dir
        return step

def test_backtracking_line_search():
    geom = AnaPot.get_geom((0.687, 1.57, 0.))

    opt = SteepestDescent(geom, thresh="gau_tight")
    opt.run()

    cs = np.array(opt.coords)

    calc = geom.calculator
    calc.plot()#show=True)

    ax = calc.ax
    ax.plot(*cs.T[:2])

    plt.show()
