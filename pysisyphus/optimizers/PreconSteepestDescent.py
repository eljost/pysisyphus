import numpy as np

from pysisyphus.line_searches import Backtracking, HagerZhang, StrongWolfe
from pysisyphus.optimizers.Optimizer import Optimizer


class PreconSteepestDescent(Optimizer):

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
        self.line_search_cls = ls_cls[self.line_search]

        self.alpha_prev = None

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        step_dir = forces / np.linalg.norm(forces)

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
