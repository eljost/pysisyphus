#!/usr/bin/env python3

from abc import abstractmethod

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.hessian_updates import bfgs_update, flowchart_update


class HessianOptimizer(Optimizer):
    hessian_update_funcs = {
        "bfgs": bfgs_update,
        "flowchart": flowchart_update,
    }

    def __init__(self, geometry, trust_radius=0.3, trust_update=True,
                 trust_min=0.01, trust_max=1, hessian_update="bfgs",
                 hessian_init="init", **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius
        self.trust_update = trust_update
        self.trust_min = trust_min
        self.trust_max = trust_max
        self.hessian_update = hessian_update
        self.hessian_update_func = self.hessian_update_funcs[hessian_update]
        self.hessian_init = hessian_init

        self.predicted_energy_changes = list()
        self.actual_energy_changes = list()

    def prepare_opt(self):
        hess_funcs = {
            # Calculate true hessian
            "calc": (self.geometry.hessian, "calculated exact"),
            # Approximate hessian
            "guess": (self.geometry.get_initial_hessian(), "approximate guess"),
            # Unit hessian
            "unit": (np.eye(self.geometry.coords.size), "unit")
        }
        try:
            self.H, hess_str = hess_funcs[self.hessian_init]
            self.log(f"Using {hess_str} hessian.")
        except KeyError:
            self.log(f"Trying to load saved hessian from '{self.hessian_init}'.")
            self.H = np.loadtxt(self.hessian_init)

    def update_trust_radius(self, coeff, last_step_norm):
        # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius/4,
                                    self.trust_min)
            self.log("Decreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
            self.trust_radius = min(self.trust_radius*2,
                                    self.trust_max)
            self.log("Increasing trust radius.")
        else:
            self.log(f"Keeping current trust radius at {self.trust_radius:.6f}")
            return
        self.log(f"Updated trust radius: {self.trust_radius:.6f}")

    def update_hessian(self):
        dx = self.steps[-1]
        dg = -(self.forces[-1] - self.forces[-2])
        dH, key = self.hessian_update_func(self.H, dx, dg)
        self.H = self.H + dH
        self.log(f"Did {key} hessian update.")

    @abstractmethod
    def optimize(self):
        pass
