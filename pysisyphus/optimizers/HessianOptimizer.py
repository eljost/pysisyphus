#!/usr/bin/env python3

from abc import abstractmethod

import numpy as np

from pysisyphus.optimizers.guess_hessians import (fischer_guess,
                                                  lindh_guess,
                                                  simple_guess,
                                                  swart_guess,
                                                  xtb_hessian,)
from pysisyphus.optimizers.hessian_updates import (bfgs_update,
                                                   flowchart_update,
                                                   damped_bfgs_update,
                                                   multi_step_update,
                                                   bofill_update,)
from pysisyphus.optimizers import line_search2
from pysisyphus.optimizers.line_search2 import poly_line_search
from pysisyphus.optimizers.Optimizer import Optimizer


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
                 hessian_recalc=None, hessian_xtb=False,
                 small_eigval_thresh=1e-8, line_search=False, **kwargs):
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
        self.hessian_xtb = hessian_xtb
        self.small_eigval_thresh = float(small_eigval_thresh)
        self.line_search = line_search
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
            # We allow only loading of cartesian hessians
            self.log(f"Trying to load saved hessian from '{self.hessian_init}'.")
            self.geometry.cart_hessian = np.loadtxt(self.hessian_init)
            # Use the previously set hessian in whatever coordinate system we
            # actually employ.
            self.H = self.geometry.hessian
        if self.hessian_init == "calc":
            hess_fn = self.get_path_for_fn("calculated_init_cart_hessian")
            # Save the cartesian hessian, as it is independent of the
            # actual coordinate system that is used.
            np.savetxt(hess_fn, self.geometry._hessian)
            self.log(f"Wrote calculated cartesian hessian to '{hess_fn}'")

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
        # Only report an unexpected increase if we actually predicted a
        # decrease.
        if (actual_change > 0) and (predicted_change < 0):
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

        # If actual and predicted energy change have different signs
        # coeff will be negative and lead to a decreased trust radius,
        # which is fine.
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius/4,
                                    self.trust_min)
            self.log("Decreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        # elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
        #
        # Only increase trust radius if last step norm corresponded approximately
        # to the trust radius.
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
            # Use xtb hessian
            self.log("Requested hessian recalculation.")
            if self.hessian_xtb:
                self.H = xtb_hessian(self.geometry)
                key = "xtb"
            # Calculated hessian at actual level of theory
            else:
                self.H = self.geometry.hessian
                key = "exact"
            if not (self.cur_cycle == 0):
                self.log(f"Recalculated {key} hessian in cycle {self.cur_cycle}.")
        elif self.hessian_multi_update:
            gradients = -np.array(self.forces)
            self.H = multi_step_update(self.H, self.steps, gradients, self.energies)
        else:
            dx = self.steps[-1]
            dg = -(self.forces[-1] - self.forces[-2])
            dH, key = self.hessian_update_func(self.H, dx, dg)
            self.H = self.H + dH
            self.log(f"Did {key} hessian update.")

    def poly_line_search(self):
        # Current energy & gradient are already appended.
        cur_energy = self.energies[-1]
        prev_energy = self.energies[-2]
        energy_increased = (cur_energy - prev_energy) > 0.

        prev_step = self.steps[-1]
        cur_grad = -self.forces[-1]
        prev_grad = -self.forces[-2]

        # TODO: always call line_search? Right now we could get some acceleration
        # as we also accept steps > 1.
        # if not energy_increased:
            # return cur_grad

        # Generate directional gradients by projecting them on the previous step.
        prev_grad_proj = prev_step @ prev_grad
        cur_grad_proj =  prev_step @ cur_grad
        cubic_result = line_search2.cubic_fit(prev_energy, cur_energy,
                                              prev_grad_proj, cur_grad_proj)
        quartic_result = line_search2.quartic_fit(prev_energy, cur_energy,
                                              prev_grad_proj, cur_grad_proj)
        # TODO: add quintic

        prev_coords = self.coords[-2]
        cur_coords = self.coords[-1]
        accept = {
            "cubic": lambda x: (x > 2.) and (x < 1),
            "quartic": lambda x: (x > 0.) and (x <= 2),
        }
        fit_result = None
        if quartic_result and accept["quartic"](quartic_result.x):
            fit_result = quartic_result
            deg = "quartic"
        elif cubic_result and accept["cubic"](cubic_result.x):
            fit_result = cubic_result
            deg = "cubic"
        # else:
            # Midpoint fallback as described by gaussian?

        if fit_result and fit_result.y < prev_energy:

            x = fit_result.x
            y = fit_result.y
            self.log(f"Did {deg} interpolation with x={x:.6f}.")
            fit_step = x * prev_step
            # Interpolate coordinates and gradient
            fit_coords = prev_coords + fit_step
            fit_grad = (1-x)*prev_grad + x*cur_grad

            # TODO: update step and other saved entries?!
            self.geometry.coords = fit_coords
            self.forces[-1] = -fit_grad
            self.energies[-1] = y
            self.coords[-1] = fit_coords.copy()
            self.cart_coords[-1] = self.geometry.cart_coords.copy()
            self.steps[-1] = fit_step
            cur_grad = fit_grad

            # self.update_hessian()
        return cur_grad

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
        self.log(f"\tnu_{verbose}={nu:.4e}")
        # Scale eigenvector so that its last element equals 1. The
        # final is step is the scaled eigenvector without the last element.
        step = step_nu[:-1] / nu
        eigval = eigenvalues[ind]
        self.log(f"\teigenvalue_{verbose}={eigval:.4e}")
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

    def log_negative_eigenvalues(self, eigvals, pre_str=""):
        neg_inds = eigvals < -self.small_eigval_thresh
        neg_eigval_str = np.array2string(eigvals[neg_inds], precision=6)
        self.log(f"{pre_str}hessian has {neg_inds.sum()} negative eigenvalue(s).")
        self.log(f"\t{neg_eigval_str}")

    @abstractmethod
    def optimize(self):
        pass
