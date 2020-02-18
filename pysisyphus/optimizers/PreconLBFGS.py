#!/usr/bin/env python3

from collections import deque

import numpy as np
from scipy.sparse.linalg import spsolve

from pysisyphus.line_searches import Backtracking, Dummy, HagerZhang, StrongWolfe
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.precon import precon_getter


class PreconLBFGS(Optimizer):

    def __init__(self, geometry, alpha_init=1., history=7, precon=True,
                 line_search="armijo", **kwargs):
        assert geometry.coord_type == "cart", \
            "Preconditioning makes only sense with 'coord_type: cart'"
        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init
        self.history = history
        self.precon = precon
        self.line_search = line_search

        ls_cls = {
            "armijo": Backtracking,
            "strong_wolfe": StrongWolfe,
            "hz": HagerZhang,
            "dummy": Dummy,
        }
        self.line_search_cls = ls_cls[self.line_search]

        self.grad_diffs = deque(maxlen=self.history)
        self.steps_ = deque(maxlen=self.history)

    def prepare_opt(self):
        if self.precon:
            self.precon_getter = precon_getter(self.geometry)

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        # Steepest descent fallback
        step = forces

        # Construct preconditoner if requested
        P = None
        if self.precon:
            P = self.precon_getter(self.geometry.coords)
            step = spsolve(P, forces)

        if self.cur_cycle > 0:
            self.grad_diffs.append(-forces - -self.forces[-2])
            self.steps_.append(self.steps[-1])
            step = -bfgs_multiply(self.steps_, self.grad_diffs, forces, P=P)

        step_dir = step / np.linalg.norm(step)

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
        return step
