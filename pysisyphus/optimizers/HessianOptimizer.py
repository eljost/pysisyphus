#!/usr/bin/env python3

from abc import abstractmethod

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.hessian_updates import (bfgs_update, flowchart_update,
                                                   damped_bfgs_update,
                                                   multi_step_update,
                                                   bofill_update)
from pysisyphus.optimizers.guess_hessians import (fischer_guess, lindh_guess,
                                                  simple_guess, swart_guess,
                                                  xtb_hessian,)


class HessianOptimizer(Optimizer):
    hessian_update_funcs = {
        "bfgs": bfgs_update,
        "damped_bfgs": damped_bfgs_update,
        "flowchart": flowchart_update,
        "mult": multi_step_update,
        "bofill": bofill_update,
    }

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(self, geometry, trust_radius=0.5, trust_update=True,
                 trust_min=0.1, trust_max=1, hessian_update="bfgs",
                 hessian_multi_update=False, hessian_init="fischer",
                 hessian_recalc=None, small_eigval_thresh=1e-8, **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_update = bool(trust_update)
        assert trust_min <= trust_max, \
                "trust_min must be <= trust_max!"
        self.trust_min = float(trust_min)
        self.trust_max = float(trust_max)
        # Constrain initial trust radius if trust_max > trust_radius
        self.trust_radius = min(trust_radius, trust_max)
        self.hessian_update = hessian_update
        self.hessian_update_func = self.hessian_update_funcs[hessian_update]
        self.hessian_multi_update = hessian_multi_update
        if self.hessian_multi_update:
            raise Exception("hessian_multi_update=True doesn't work yet!")
        self.hessian_init = hessian_init
        self.hessian_recalc = hessian_recalc
        self.small_eigval_thresh = float(small_eigval_thresh)
        assert self.small_eigval_thresh > 0., "small_eigval_thresh must be > 0.!"

        # Allow only calculated or unit hessian for geometries that don't
        # use internal coordinates.
        if (not hasattr(self.geometry, "internal")
            or (self.geometry.internal is None)):
            if self.hessian_init != "calc":
                self.hessian_init = "unit"

        self.predicted_energy_changes = list()

    def prepare_opt(self):
        # We use lambdas to avoid premature evaluation of the dict items.
        hess_funcs = {
            # Calculate true hessian
            "calc": lambda: (self.geometry.hessian, "calculated exact"),
            # Unit hessian
            "unit": lambda: (np.eye(self.geometry.coords.size), "unit"),
            # Fischer model hessian
            "fischer": lambda: (fischer_guess(self.geometry), "Fischer"),
            # Lindh model hessian
            "lindh": lambda: (lindh_guess(self.geometry), "Lindh"),
            # Simple (0.5, 0.2, 0.1) model hessian
            "simple": lambda: (simple_guess(self.geometry), "simple"),
            # Swart model hessian
            "swart": lambda: (swart_guess(self.geometry), "Swart"),
            # XTB hessian
            "xtb": lambda: (xtb_hessian(self.geometry), "XTB"),
        }
        try:
            self.H, hess_str = hess_funcs[self.hessian_init]()
            self.log(f"Using {hess_str} hessian.")
        except KeyError:
            self.log(f"Trying to load saved hessian from '{self.hessian_init}'.")
            self.H = np.loadtxt(self.hessian_init)
        if self.hessian_init == "calc":
            hess_fn = self.get_path_for_fn("calculated_init_hessian")
            np.savetxt(hess_fn, self.H)
            self.log(f"Wrote calculated hessian to '{hess_fn}'")

        if (hasattr(self.geometry, "coord_type")
            and self.geometry.coord_type == "dlc"):
            U = self.geometry.internal.U
            self.H = U.T.dot(self.H).dot(U)

    def update_trust_radius(self):
        # The predicted change should be calculated at the end of optimize
        # of the previous cycle.
        assert len(self.predicted_energy_changes) == len(self.forces)-1, \
            "Did you forget to append to self.predicted_energy_changes?"
        predicted_change = self.predicted_energy_changes[-1]
        actual_change = self.energies[-1] - self.energies[-2]
        if actual_change > 0:
            print(f"Energy increased by {actual_change:.6f} au! " \
                  f"Cur. trust={self.trust_radius:.6f}.")
            self.log(f"Energy increased by {actual_change:.6f} au!")
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
        # elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
        elif coeff > 0.75 and abs(self.trust_radius - last_step_norm) <= 1e-3:
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
        elif self.hessian_multi_update:
            gradients = -np.array(self.forces)
            self.H = multi_step_update(self.H, self.steps, gradients, self.energies)
        else:
            dx = self.steps[-1]
            dg = -(self.forces[-1] - self.forces[-2])
            dH, key = self.hessian_update_func(self.H, dx, dg)
            self.H = self.H + dH
            self.log(f"Did {key} hessian update.")

    def solve_rfo(self, rfo_mat, kind="min"):
        eigenvalues, eigenvectors = np.linalg.eig(rfo_mat)
        eigenvalues = eigenvalues.real
        eigenvectors = eigenvectors.real
        sorted_inds = np.argsort(eigenvalues)

        # Depending on wether we want to minimize (maximize) along
        # the mode(s) in the rfo mat we have to select the smallest
        # (biggest) eigenvalue and corresponding eigenvector.
        first_or_last, verbose = self.rfo_dict[kind]
        ind = sorted_inds[first_or_last]
        # Given sorted eigenvalue-indices (sorted_inds) use the first
        # (smallest eigenvalue) or the last (largest eigenvalue) index.
        step_nu = eigenvectors.T[ind]
        nu = step_nu[-1]
        self.log(f"nu_{verbose}={nu:.4e}")
        # Scale eigenvector so that its last element equals 1. The
        # final is step is the scaled eigenvector without the last element.
        step = step_nu[:-1] / nu
        eigval = eigenvalues[ind]
        self.log(f"eigenvalue_{verbose}={eigval:.4e}")
        return step, eigval, nu

    def filter_small_eigvals(self, eigvals, eigvecs):
        small_inds = np.abs(eigvals) < self.small_eigval_thresh
        eigvals = eigvals[~small_inds]
        eigvecs = eigvecs[:,~small_inds]
        small_num = sum(small_inds)
        self.log(f"Found {small_num} small eigenvalues in hessian. Removed "
                  "corresponding eigenvalues and eigenvectors.")
        assert small_num <= 6, \
             "Expected at most 6 small eigenvalues in cartesian hessian " \
            f"but found {small_num}!"
        return eigvals, eigvecs

    @abstractmethod
    def optimize(self):
        pass
