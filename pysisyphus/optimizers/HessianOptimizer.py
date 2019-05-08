#!/usr/bin/env python3

from abc import abstractmethod

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.hessian_updates import bfgs_update, flowchart_update
from pysisyphus.optimizers.guess_hessians import lindh_guess


class HessianOptimizer(Optimizer):
    hessian_update_funcs = {
        "bfgs": bfgs_update,
        "flowchart": flowchart_update,
    }

    def __init__(self, geometry, trust_radius=0.5, trust_update=True,
                 trust_min=0.01, trust_max=1, hessian_update="bfgs",
                 hessian_init="guess", hessian_recalc=None, **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius
        self.trust_update = trust_update
        self.trust_min = trust_min
        self.trust_max = trust_max
        self.hessian_update = hessian_update
        self.hessian_update_func = self.hessian_update_funcs[hessian_update]
        self.hessian_init = hessian_init
        self.hessian_recalc = hessian_recalc

        self.predicted_energy_changes = list()

    def prepare_opt(self):
        # We use lambdas to avoid premature evaluation of the dict items.
        hess_funcs = {
            # Calculate true hessian
            "calc": lambda: (self.geometry.hessian, "calculated exact"),
            # Approximate hessian
            "guess": lambda: (self.geometry.get_initial_hessian(), "approximate guess"),
            # Unit hessian
            "unit": lambda: (np.eye(self.geometry.coords.size), "unit"),
            # Lindh model hessian
            "lindh": lambda: (lindh_guess(self.geometry), "Lindh"),
        }
        try:
            self.H, hess_str = hess_funcs[self.hessian_init]()
            self.log(f"Using {hess_str} hessian.")
        except KeyError:
            self.log(f"Trying to load saved hessian from '{self.hessian_init}'.")
            self.H = np.loadtxt(self.hessian_init)
        if self.hessian_init == "calc":
            hess_fn = "calculated_init_hessian"
            np.savetxt(hess_fn, self.H)
            self.log(f"Wrote calculated hessian to '{hess_fn}'")

    def update_trust_radius(self):
        assert len(self.predicted_energy_changes) == len(self.forces)-1, \
            "Did you forget to append to self.predicted_energy_changes?"
        predicted_change = self.predicted_energy_changes[-1]
        actual_change = self.energies[-1] - self.energies[-2]
        coeff = actual_change / predicted_change
        self.log(f"Predicted change: {predicted_change:.4e} au")
        self.log(f"Actual change: {actual_change:.4e} au")
        self.log(f"Coefficient: {coeff:.2%}")
        if self.trust_update:
            step = self.steps[-1]
            last_step_norm = np.linalg.norm(step)
            self.get_new_trust_radius(coeff, last_step_norm)
        else:
            self.log("Skipping trust radius update")

    def get_new_trust_radius(self, coeff, last_step_norm):
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
        if self.hessian_recalc and (self.cur_cycle % self.hessian_recalc) == 0:
            self.H = self.geometry.hessian
            if not (self.cur_cycle == 0):
                self.log(f"Recalculated exact hessian in cycle {self.cur_cycle}.")
        else:
            dx = self.steps[-1]
            dg = -(self.forces[-1] - self.forces[-2])
            dH, key = self.hessian_update_func(self.H, dx, dg)
            self.H = self.H + dH
            self.log(f"Did {key} hessian update.")

    @abstractmethod
    def optimize(self):
        pass
