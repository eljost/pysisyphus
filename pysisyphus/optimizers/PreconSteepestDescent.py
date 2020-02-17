import numpy as np
from scipy.sparse.linalg import spsolve

from pysisyphus.line_searches import Backtracking, HagerZhang, StrongWolfe
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.precon import precon_getter


class PreconSteepestDescent(Optimizer):

    def __init__(self, geometry, alpha_init=0.5, line_search="armijo",
                 precon=True, **kwargs):
        assert geometry.coord_type == "cart", \
            "Preconditioning makes only sense with 'coord_type: cart'"
        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init
        self.line_search = line_search
        self.precon = precon

        ls_cls = {
            "armijo": Backtracking,
            "strong_wolfe": StrongWolfe,
            "hz": HagerZhang,
        }
        self.line_search_cls = ls_cls[self.line_search]

        self.alpha_prev = None

    def prepare_opt(self):
        if self.precon:
            self.precon_getter = precon_getter(self.geometry)

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        # Plain steepst descent (fallback)
        step = forces

        # Preconditoned steepest descent if requested
        if self.precon:
            P = self.precon_getter(self.geometry.coords)
            # Solve:
            #   step = -P^-1 * grad
            #   P step = -grad
            step = spsolve(P, forces)

        step_dir = step / np.linalg.norm(step)

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
