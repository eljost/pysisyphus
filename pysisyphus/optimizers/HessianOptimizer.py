#!/usr/bin/env python3

import numpy as np

from pysisyphus.intcoords.helpers import get_step
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
                 hessian_recalc=None, hessian_recalc_adapt=None, hessian_xtb=False,
                 small_eigval_thresh=1e-8, line_search=False,
                 alpha0=1., max_micro_cycles=25,
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_update = bool(trust_update)
        assert trust_min <= trust_max, \
                "trust_min must be <= trust_max!"
        self.trust_min = float(trust_min)
        self.trust_max = float(trust_max)
        # Constrain initial trust radius if trust_max > trust_radius
        self.trust_radius = min(trust_radius, trust_max)
        self.log(f"Initial trust radius: {self.trust_radius:.6f}")
        self.hessian_update = hessian_update
        self.hessian_update_func = self.hessian_update_funcs[hessian_update]
        self.hessian_multi_update = hessian_multi_update
        if self.hessian_multi_update:
            raise Exception("hessian_multi_update=True doesn't work yet!")
        self.hessian_init = hessian_init
        self.hessian_recalc = hessian_recalc
        self.hessian_recalc_adapt = hessian_recalc_adapt
        self.hessian_xtb = hessian_xtb
        self.small_eigval_thresh = float(small_eigval_thresh)
        self.line_search = bool(line_search)
        # Restricted-step related
        self.alpha0 = alpha0
        self.max_micro_cycles = int(max_micro_cycles)
        assert max_micro_cycles >= 1

        assert self.small_eigval_thresh > 0., "small_eigval_thresh must be > 0.!"
        self.hessian_recalc_in = None
        self.adapt_norm = None

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

        if self.hessian_recalc_adapt:
            self.adapt_norm = np.linalg.norm(self.geometry.forces)

        if self.hessian_recalc:
            # Already substract one, as we don't do a hessian update in
            # the first cycle.
            self.hessian_recalc_in = self.hessian_recalc - 1

    def update_trust_radius(self):
        # The predicted change should be calculated at the end of optimize
        # of the previous cycle.
        assert len(self.predicted_energy_changes) == len(self.forces)-1, \
            "Did you forget to append to self.predicted_energy_changes?"
        self.log("Trust radius update")
        self.log(f"\tCurrent trust radius: {self.trust_radius:.6f}")
        predicted_change = self.predicted_energy_changes[-1]
        actual_change = self.energies[-1] - self.energies[-2]
        # Only report an unexpected increase if we actually predicted a
        # decrease.
        if (actual_change > 0) and (predicted_change < 0):
            print(f"Energy increased by {actual_change:.6f} au! " \
                  f"Cur. trust={self.trust_radius:.6f}.")
            self.log(f"Energy increased by {actual_change:.6f} au!")
        coeff = actual_change / predicted_change
        self.log(f"\tPredicted change: {predicted_change:.4e} au")
        self.log(f"\tActual change: {actual_change:.4e} au")
        self.log(f"\tCoefficient: {coeff:.2%}")
        if self.trust_update:
            step = self.steps[-1]
            last_step_norm = np.linalg.norm(step)
            self.set_new_trust_radius(coeff, last_step_norm)
        else:
            self.log("\tSkipped trust radius update")

    def set_new_trust_radius(self, coeff, last_step_norm):
        # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1

        # If actual and predicted energy change have different signs
        # coeff will be negative and lead to a decreased trust radius,
        # which is fine.
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius/4,
                                    self.trust_min)
            self.log("\tDecreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        # elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
        #
        # Only increase trust radius if last step norm corresponded approximately
        # to the trust radius.
        elif coeff > 0.75 and abs(self.trust_radius - last_step_norm) <= 1e-3:
            self.trust_radius = min(self.trust_radius*2,
                                    self.trust_max)
            self.log("\tIncreasing trust radius.")
        else:
            self.log(f"\tKeeping current trust radius at {self.trust_radius:.6f}")
            return
        self.log(f"\tUpdated trust radius: {self.trust_radius:.6f}")

    def update_hessian(self):
        # Compare current forces to reference forces to see if we shall recalc the
        # hessian.
        try:
            cur_norm = np.linalg.norm(self.forces[-1])
            ref_norm = self.adapt_norm / self.hessian_recalc_adapt
            recalc_adapt = cur_norm <= ref_norm
            self.log( "Check for adaptive hessian recalculation: "
                     f"{cur_norm:.6f} <= {ref_norm:.6f}, {recalc_adapt}"
            )
        except TypeError:
            recalc_adapt = False

        try:
            self.hessian_recalc_in = max(self.hessian_recalc_in-1, 0)
            self.log(f"Recalculation of hessian in {self.hessian_recalc_in} cycle(s).")
        except TypeError:
            self.hessian_recalc_in = None

        # Update reference norm if needed
        # TODO: Decide on whether to update the norm when the recalculation is
        # initiated by 'recalc'.
        if recalc_adapt:
            self.adapt_norm = cur_norm

        recalc = (self.hessian_recalc_in == 0)

        if recalc or recalc_adapt:
            # Use xtb hessian
            self.log("Requested hessian recalculation.")
            if self.hessian_xtb:
                self.H = xtb_hessian(self.geometry)
                key = "xtb"
            # Calculated hessian at actual level of theory
            else:
                self.H = self.geometry.hessian
                key = "exact"
                hess_fn = self.get_path_for_fn(
                            f"calculated_cart_hessian_cycle_{self.cur_cycle}"
                )
                np.savetxt(hess_fn, self.geometry._hessian)
                self.log(f"Wrote calculated cartesian hessian to '{hess_fn}'")
            if not (self.cur_cycle == 0):
                self.log(f"Recalculated {key} hessian in cycle {self.cur_cycle}.")
            # Reset counter. It is also reset when the recalculation was initiated
            # by the adaptive formulation.
            self.hessian_recalc_in = self.hessian_recalc
        # elif self.hessian_multi_update:
            # gradients = -np.array(self.forces)
            # self.H = multi_step_update(self.H, self.steps, gradients, self.energies)
        # Simple hessian update
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
        # energy_increased = (cur_energy - prev_energy) > 0.

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
        prev_coords = self.coords[-2]
        accept = {
            # cubic is disabled for now as it does not seem to help
            "cubic": lambda x: (x > 2.) and (x < 1),  # lgtm [py/redundant-comparison]
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

        fit_energy = None
        fit_grad = None
        fit_coords = None
        fit_step = None
        if fit_result and fit_result.y < prev_energy:
            x = fit_result.x
            fit_energy = fit_result.y
            self.log(f"Did {deg} interpolation with x={x:.6f}.")

            # Interpolate coordinates and gradient
            fit_step = x * prev_step
            fit_coords = prev_coords + fit_step
            # The commented lines below would be correct if we would want
            # the step from the previous coordinates and not the current ones.
            # fit_step = (1-x) * -prev_step
            # fit_coords = cur_coords + fit_step
            fit_grad = (1-x)*prev_grad + x*cur_grad
        return fit_energy, fit_grad, fit_coords, fit_step

    def poly_line_search_v2(self, hessian=None):
        assert len(self.energies) == len(self.coords) == len(self.forces)
        # Find previous best point
        prev_best_ind = np.argmin(self.energies[:-1])
        prev_best_energy = self.energies[prev_best_ind]
        prev_best_coords = self.coords[prev_best_ind]
        prev_best_grad = -self.forces[prev_best_ind]

        # Current point. Current energy & gradient are already appended.
        cur_energy = self.energies[-1]
        cur_grad = -self.forces[-1]

        at_best_energy = cur_energy < prev_best_energy

        # This wont work for internals
        # step = prev_coords - cur_coords
        # This should work for all coord_types
        # tmp_geom = self.geometry.copy(check_bends=False)
        # tmp_geom.coords = prev_best_coords
        # step = tmp_geom - self.geometry
        step = get_step(self.geometry, prev_best_coords)

        if hessian is not None:
            hess_proj = float(step[None,:].dot(hessian).dot(step[:,None]))

        # Generate directional gradients by projecting them on the previous step.
        prev_grad_proj = step @ prev_best_grad
        cur_grad_proj =  step @ cur_grad
        cubic_result = line_search2.cubic_fit(prev_best_energy, cur_energy,
                                              prev_grad_proj, cur_grad_proj)
        quartic_result = line_search2.quartic_fit(prev_best_energy, cur_energy,
                                                  prev_grad_proj, cur_grad_proj)
        quintic_result = None
        if hessian is not None:
            quintic_result = line_search2.quintic_fit(prev_best_energy, cur_energy,
                                                      prev_grad_proj, cur_grad_proj,
                                                      hess_proj, hess_proj)
        accept_if_best = lambda x: True if at_best_energy else (0. < x < 1.)
        accept = {
            "cubic": lambda x: 0. < x < 1.,
            # "quartic": accept_if_best,
            "quartic": lambda x: 0. < x < 2.,
            "quintic": accept_if_best,
        }
        fit_result = None
        # import pdb; pdb.set_trace()
        if quintic_result and accept["quintic"](quintic_result.x):
            fit_result = quintic_result
            deg = "quintic"
        elif quartic_result and accept["quartic"](quartic_result.x):
            fit_result = quartic_result
            deg = "quartic"
        elif cubic_result and accept["cubic"](cubic_result.x):
            fit_result = cubic_result
            deg = "cubic"
        # else:
            # Midpoint fallback as described by gaussian?

        fit_energy = None
        fit_grad = None
        fit_coords = None
        fit_step = None
        if fit_result and fit_result.y < prev_best_energy:
            x = fit_result.x
            fit_energy = fit_result.y
            self.log(f"Did {deg} interpolation with x={x:.6f}.")

            # Interpolate coordinates and gradient
            fit_step = x * step
            fit_coords = prev_best_coords + fit_step
            # The commented lines below would be correct if we would want
            # the step from the previous coordinates and not the current ones.
            # fit_step = (1-x) * -prev_step
            # fit_coords = cur_coords + fit_step
            fit_grad = (1-x)*prev_best_grad + x*cur_grad
        return fit_energy, fit_grad, fit_coords, fit_step

    def solve_rfo(self, rfo_mat, kind="min"):
        self.log("Diagonalizing augmented Hessian:")
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
        # TODO: Root following like in optking?!
        nu = step_nu[-1]
        self.log(f"\tnu_{verbose}={nu:.4e}")
        # Scale eigenvector so that its last element equals 1. The
        # final is step is the scaled eigenvector without the last element.
        step = step_nu[:-1] / nu
        eigval = eigenvalues[ind]
        self.log(f"\teigenvalue_{verbose}={eigval:.4e}")
        return step, eigval, nu

    def filter_small_eigvals(self, eigvals, eigvecs, mask=False):
        small_inds = np.abs(eigvals) < self.small_eigval_thresh
        eigvals = eigvals[~small_inds]
        eigvecs = eigvecs[:,~small_inds]
        small_num = sum(small_inds)
        self.log(f"Found {small_num} small eigenvalues in hessian. Removed "
                  "corresponding eigenvalues and eigenvectors.")
        assert small_num <= 6, \
             "Expected at most 6 small eigenvalues in cartesian hessian " \
            f"but found {small_num}!"
        if mask:
            return eigvals, eigvecs, small_inds
        else:
            return eigvals, eigvecs

    def log_negative_eigenvalues(self, eigvals, pre_str=""):
        neg_inds = eigvals < -self.small_eigval_thresh
        neg_eigval_str = np.array2string(eigvals[neg_inds], precision=6)
        self.log(f"{pre_str}hessian has {neg_inds.sum()} negative eigenvalue(s).")
        self.log(f"\t{neg_eigval_str}")

    def housekeeping(self):
        """Calculate gradient and energy. Update trust radius and hessian
        if needed. Return energy, gradient and hessian for the current cycle."""
        gradient = self.geometry.gradient
        energy = self.geometry.energy
        self.forces.append(-gradient)
        self.energies.append(energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

        H = self.H
        if self.geometry.internal:
            # Shift eigenvalues of orthogonal part to high values, so they
            # don't contribute to the actual step.
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection may break it?!
            H = (H_proj + H_proj.T) / 2

        eigvals, eigvecs = np.linalg.eigh(H)
        # Neglect small eigenvalues
        eigvals, eigvecs = self.filter_small_eigvals(eigvals, eigvecs)

        return energy, gradient, H, eigvals, eigvecs

    def get_augmented_hessian(self, eigvals, gradient, alpha=1.):
        dim_ = eigvals.size + 1
        H_aug = np.zeros((dim_, dim_))
        H_aug[:dim_-1,:dim_-1] = np.diag(eigvals/alpha)
        H_aug[-1,:-1] = gradient
        H_aug[:-1,-1] = gradient

        H_aug[:-1,-1] /= alpha

        return H_aug

    def get_alpha_step(self, cur_alpha, rfo_eigval, step_norm, eigvals, gradient):
        # Derivative of the squared step w.r.t. alpha
        numer = gradient**2
        denom = (eigvals - rfo_eigval * cur_alpha)**3
        quot = np.sum(numer / denom)
        self.log(f"quot={quot:.6f}")
        dstep2_dalpha = (2*rfo_eigval/(1+step_norm**2 * cur_alpha)
                         * np.sum(gradient**2
                                  / ((eigvals - rfo_eigval * cur_alpha)**3)
                           )
        )
        self.log(f"analytic deriv.={dstep2_dalpha:.6f}")
        # Update alpha
        alpha_step = (2*(self.trust_radius*step_norm - step_norm**2)
                      / dstep2_dalpha
        )
        self.log(f"alpha_step={alpha_step:.4f}")
        assert (cur_alpha + alpha_step) > 0, "alpha must not be negative!"
        return alpha_step

    def get_rs_step(self, eigvals, eigvecs, gradient, name="RS"):
        # Transform gradient to basis of eigenvectors
        gradient_ = eigvecs.T.dot(gradient)

        alpha = self.alpha0
        for mu in range(self.max_micro_cycles):
            self.log(f"{name} micro cycle {mu:02d}, alpha={alpha:.6f}")
            H_aug = self.get_augmented_hessian(eigvals, gradient_, alpha)
            rfo_step_, eigval_min, nu = self.solve_rfo(H_aug, "min")
            rfo_norm_ = np.linalg.norm(rfo_step_)
            self.log(f"norm(rfo step)={rfo_norm_:.6f}")

            if (rfo_norm_ < self.trust_radius) or abs(rfo_norm_ - self.trust_radius) <= 1e-3:
                step_ = rfo_step_
                break

            alpha_step = self.get_alpha_step(alpha, eigval_min, rfo_norm_, eigvals, gradient_)
            alpha += alpha_step
            self.log("")
        else:
            self.log( "RS algorithm did not produce a desired step length "
                     f"after {self.max_micro_cycles} micro cycles. Using "
                      "simple downscaled step with alpha=1."
            )
            H_aug = self.get_augmented_hessian(eigvals, gradient_, alpha=1.)
            rfo_step_, eigval_min, nu = self.solve_rfo(H_aug, "min")
            rfo_norm_ = np.linalg.norm(rfo_step_)
            # This should always be True if the above algorithm failed but we
            # keep this line here nonetheless to make it more obvious what
            # we are doing.
            if rfo_norm_ > self.trust_radius:
                step_ = rfo_step_ / rfo_norm_ * self.trust_radius

        # Transform step back to original basis
        step = eigvecs.dot(step_)
        return step
