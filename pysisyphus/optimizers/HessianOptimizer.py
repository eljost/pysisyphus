from math import sqrt
from typing import Literal, Optional

import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import rms
from pysisyphus.io.hessian import save_hessian
from pysisyphus.optimizers.guess_hessians import (
    get_guess_hessian,
    xtb_hessian,
    HessInit,
)
from pysisyphus.optimizers.hessian_updates import (
    bfgs_update,
    flowchart_update,
    damped_bfgs_update,
    bofill_update,
    ts_bfgs_update,
    ts_bfgs_update_org,
    ts_bfgs_update_revised,
)
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.exceptions import OptimizationError


def dummy_hessian_update(H, dx, dg):
    return np.zeros_like(H), "no"


HESS_UPDATE_FUNCS = {
    "none": dummy_hessian_update,
    None: dummy_hessian_update,
    False: dummy_hessian_update,
    "bfgs": bfgs_update,
    "damped_bfgs": damped_bfgs_update,
    "flowchart": flowchart_update,
    "bofill": bofill_update,
    "ts_bfgs": ts_bfgs_update,
    "ts_bfgs_org": ts_bfgs_update_org,
    "ts_bfgs_rev": ts_bfgs_update_revised,
}
HessUpdate = Literal[
    "none",
    None,
    False,
    "bfgs",
    "damped_bfgs",
    "flowchart",
    "bofill",
    "ts_bfgs",
    "ts_bfgs_org",
    "ts_bfgs_rev",
]


class HessianOptimizer(Optimizer):

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(
        self,
        geometry: Geometry,
        trust_radius: float = 0.5,
        trust_update: bool = True,
        trust_min: float = 0.1,
        trust_max: float = 1,
        max_energy_incr: Optional[float] = None,
        hessian_update: HessUpdate = "bfgs",
        hessian_init: HessInit = "fischer",
        hessian_recalc: Optional[int] = None,
        hessian_recalc_adapt: Optional[float] = None,
        hessian_xtb: bool = False,
        hessian_recalc_reset: bool = False,
        small_eigval_thresh: float = 1e-8,
        line_search: bool = False,
        alpha0: float = 1.0,
        max_micro_cycles: int = 25,
        rfo_overlaps: bool = False,
        **kwargs,
    ) -> None:
        """Baseclass for optimizers utilizing Hessian information.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        trust_radius
            Initial trust radius in whatever unit the optimization is carried out.
        trust_update
            Whether to update the trust radius throughout the optimization.
        trust_min
            Minimum trust radius.
        trust_max
            Maximum trust radius.
        max_energy_incr
            Maximum allowed energy increased after a faulty step. Optimization is
            aborted when the threshold is exceeded.
        hessian_update
            Type of Hessian update. Defaults to BFGS for minimizations and Bofill
            for saddle point searches.
        hessian_init
            Type of initial model Hessian.
        hessian_recalc
            Recalculate exact Hessian every n-th cycle instead of updating it.
        hessian_recalc_adapt
            Use a more flexible scheme to determine Hessian recalculation. Undocumented.
        hessian_xtb
            Recalculate the Hessian at the GFN2-XTB level of theory.
        hessian_recalc_reset
            Whether to skip Hessian recalculation after reset. Undocumented.
        small_eigval_thresh
            Threshold for small eigenvalues. Eigenvectors belonging to eigenvalues
            below this threshold are discardewd.
        line_search
            Whether to carry out a line search. Not implemented by a subclassing
            optimizers.
        alpha0
            Initial alpha for restricted-step (RS) procedure.
        max_micro_cycles
            Maximum number of RS iterations.
        rfo_overlaps
            Enable mode-following in RS procedure.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Optimizer baseclass.
        """
        super().__init__(geometry, **kwargs)

        assert not issubclass(
            type(geometry), ChainOfStates
        ), "HessianOptimizer can't be used for and ChainOfStates objects!"

        self.trust_update = bool(trust_update)
        assert trust_min <= trust_max, "trust_min must be <= trust_max!"
        self.trust_min = float(trust_min)
        self.trust_max = float(trust_max)
        self.max_energy_incr = max_energy_incr
        # Constrain initial trust radius if trust_max > trust_radius
        self.trust_radius = min(trust_radius, trust_max)
        self.log(f"Initial trust radius: {self.trust_radius:.6f}")
        self.hessian_update = hessian_update
        self.hessian_update_func = HESS_UPDATE_FUNCS[hessian_update]
        self.hessian_init = hessian_init
        self.hessian_recalc = hessian_recalc
        self.hessian_recalc_adapt = hessian_recalc_adapt
        self.hessian_xtb = hessian_xtb
        self.hessian_recalc_reset = hessian_recalc_reset
        self.small_eigval_thresh = float(small_eigval_thresh)
        self.line_search = bool(line_search)
        # Restricted-step related
        self.alpha0 = alpha0
        self.max_micro_cycles = int(max_micro_cycles)
        assert max_micro_cycles >= 0
        self.rfo_overlaps = rfo_overlaps

        assert self.small_eigval_thresh > 0.0, "small_eigval_thresh must be > 0.!"
        if not self.restarted:
            self.hessian_recalc_in = None
            self.adapt_norm = None
            self.predicted_energy_changes = list()
        if (
            # Allow actually calculated Hessians for all coordinate systems
            self.hessian_init not in ("calc", "xtb", "xtb1", "xtbff")
            # But disable model Hessian for Cartesian optimizations
            and self.geometry.coord_type in ("cart", "cartesian", "mwcartesian")
        ):
            self.hessian_init = "unit"
            self.log(
                f"Chosen initial (model) Hessian is incompatible with current "
                f"coord_type: {self.geometry.coord_type}!"
            )

        self._prev_eigvec_min = None
        self._prev_eigvec_max = None

    @property
    def prev_eigvec_min(self):
        return self._prev_eigvec_min

    @prev_eigvec_min.setter
    def prev_eigvec_min(self, prev_eigvec_min):
        if self.rfo_overlaps:
            self._prev_eigvec_min = prev_eigvec_min

    @property
    def prev_eigvec_max(self):
        return self._prev_eigvec_max

    @prev_eigvec_min.setter
    def prev_eigvec_max(self, prev_eigvec_max):
        if self.rfo_overlaps:
            self._prev_eigvec_max = prev_eigvec_max

    def reset(self):
        # Don't recalculate the hessian if we have to reset the optimizer
        hessian_init = self.hessian_init
        if (
            (not self.hessian_recalc_reset)
            and hessian_init == "calc"
            and self.geometry.coord_type != "cart"
        ):
            hessian_init = "fischer"
        self.prepare_opt(hessian_init)

    def save_hessian(self):
        # Don't try to save Hessians of analytical potentials
        if self.geometry.is_analytical_2d:
            return

        h5_fn = self.get_path_for_fn(f"hess_calc_cyc_{self.cur_cycle}.h5")
        # Save the cartesian hessian, as it is independent of the
        # actual coordinate system that is used.
        save_hessian(
            h5_fn,
            self.geometry,
            self.geometry.cart_hessian,
            self.geometry.energy,
            self.geometry.calculator.mult,
        )
        self.log(f"Wrote calculated cartesian Hessian to '{h5_fn}'")

    def prepare_opt(self, hessian_init=None):
        if hessian_init is None:
            hessian_init = self.hessian_init

        self.H, hess_str = get_guess_hessian(self.geometry, hessian_init)
        if self.hessian_init != "calc" and self.geometry.is_analytical_2d:
            assert self.H.shape == (3, 3)
            self.H[2, 2] = 0.0

        msg = f"Using {hess_str} Hessian"
        if hess_str == "saved":
            msg += f" from '{hessian_init}'"
        self.log(msg)

        # Dump to disk if hessian was calculated
        if self.hessian_init == "calc":
            self.save_hessian()

        if (
            hasattr(self.geometry, "coord_type")
            and self.geometry.coord_type == "dlc"
            # Calculated Hessian is already in DLC
            and hessian_init != "calc"
        ):
            U = self.geometry.internal.U
            self.H = U.T.dot(self.H).dot(U)

        if self.hessian_recalc_adapt:
            self.adapt_norm = np.linalg.norm(self.geometry.forces)

        if self.hessian_recalc:
            # Already substract one, as we don't do a hessian update in
            # the first cycle.
            self.hessian_recalc_in = self.hessian_recalc - 1

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "adapt_norm": self.adapt_norm,
            "H": self.H.tolist(),
            "hessian_recalc_in": self.hessian_recalc_in,
            "predicted_energy_changes": self.predicted_energy_changes,
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.adapt_norm = opt_restart_info["adapt_norm"]
        self.H = np.array(opt_restart_info["H"])
        self.hessian_recalc_in = opt_restart_info["hessian_recalc_in"]
        self.predicted_energy_changes = opt_restart_info["predicted_energy_changes"]

    def update_trust_radius(self):
        # The predicted change should be calculated at the end of optimize
        # of the previous cycle.
        assert (
            len(self.predicted_energy_changes) == len(self.forces) - 1
        ), "Did you forget to append to self.predicted_energy_changes?"
        self.log("Trust radius update")
        self.log(f"\tCurrent trust radius: {self.trust_radius:.6f}")
        predicted_change = self.predicted_energy_changes[-1]
        actual_change = self.energies[-1] - self.energies[-2]
        # Only report an unexpected increase if we actually predicted a
        # decrease.
        unexpected_increase = (actual_change > 0) and (predicted_change < 0)
        old_trust = self.trust_radius
        if unexpected_increase:
            self.log(f"Energy increased by {actual_change:.6f} au!")
            if self.max_energy_incr and (actual_change > self.max_energy_incr):
                raise OptimizationError("Actual energy change too high!")
        coeff = actual_change / predicted_change
        self.log(f"\tPredicted change: {predicted_change:.4e} au")
        self.log(f"\tActual change: {actual_change:.4e} au")
        self.log(f"\tCoefficient: {coeff:.2%}")
        step = self.steps[-1]
        last_step_norm = np.linalg.norm(step)
        self.set_new_trust_radius(coeff, last_step_norm)
        if unexpected_increase:
            self.table.print(
                f"Unexpected energy increase ({actual_change:.6f} au)! "
                f"Trust radius: old={old_trust:.4}, new={self.trust_radius:.4}"
            )

    def set_new_trust_radius(self, coeff, last_step_norm):
        # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1

        # If actual and predicted energy change have different signs
        # coeff will be negative and lead to a decreased trust radius,
        # which is fine.
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius / 4, self.trust_min)
            self.log("\tDecreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        # elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
        #
        # Only increase trust radius if last step norm corresponded approximately
        # to the trust radius.
        elif coeff > 0.75 and abs(self.trust_radius - last_step_norm) <= 1e-3:
            self.trust_radius = min(self.trust_radius * 2, self.trust_max)
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
            self.log(
                "Check for adaptive Hessian recalculation: "
                f"{cur_norm:.6f} <= {ref_norm:.6f}, {recalc_adapt}"
            )
        except TypeError:
            recalc_adapt = False

        try:
            self.hessian_recalc_in = max(self.hessian_recalc_in - 1, 0)
            self.log(f"Recalculation of Hessian in {self.hessian_recalc_in} cycle(s).")
        except TypeError:
            self.hessian_recalc_in = None

        # Update reference norm if needed
        # TODO: Decide on whether to update the norm when the recalculation is
        # initiated by 'recalc'.
        if recalc_adapt:
            self.adapt_norm = cur_norm

        recalc = self.hessian_recalc_in == 0

        if recalc or recalc_adapt:
            # Use xtb hessian
            self.log("Requested Hessian recalculation.")
            if self.hessian_xtb:
                self.H = xtb_hessian(self.geometry)
                key = "xtb"
            # Calculated hessian at actual level of theory
            else:
                self.H = self.geometry.hessian
                key = "exact"
                self.save_hessian()
            if not (self.cur_cycle == 0):
                self.log(f"Recalculated {key} Hessian in cycle {self.cur_cycle}.")
            # Reset counter. It is also reset when the recalculation was initiated
            # by the adaptive formulation.
            self.hessian_recalc_in = self.hessian_recalc
        # Simple hessian update
        else:
            dx = self.steps[-1]
            dg = -(self.forces[-1] - self.forces[-2])
            curv_cond = dx.dot(dg)
            if curv_cond < 0.0:
                self.log(
                    f"Curvature condition (s·y = {curv_cond:.4f} < 0) not satisfied!"
                )
            dH, key = self.hessian_update_func(self.H, dx, dg)
            self.H = self.H + dH
            self.log(f"Did {key} Hessian update.")

    def solve_rfo(self, rfo_mat, kind="min", prev_eigvec=None):
        # When using the restricted step variant of RFO the RFO matrix
        # may not be symmetric. Thats why we can't use eigh here.
        eigenvalues, eigenvectors = np.linalg.eig(rfo_mat)
        self.log("\tdiagonalized augmented Hessian")
        eigenvalues = eigenvalues.real
        eigenvectors = eigenvectors.real
        sorted_inds = np.argsort(eigenvalues)

        # Depending on wether we want to minimize (maximize) along
        # the mode(s) in the rfo mat we have to select the smallest
        # (biggest) eigenvalue and corresponding eigenvector.
        first_or_last, verbose = self.rfo_dict[kind]
        # Given sorted eigenvalue-indices (sorted_inds) use the first
        # (smallest eigenvalue) or the last (largest eigenvalue) index.
        if prev_eigvec is None:
            ind = sorted_inds[first_or_last]
        else:
            ovlps = np.array([prev_eigvec.dot(ev) for ev in eigenvectors.T])
            naive_ind = sorted_inds[first_or_last]
            ind = np.abs(ovlps).argmax()
            self.log(
                f"Overlap: {ind} ({eigenvalues[ind]:.6f}), "
                f"Naive: {naive_ind} ({eigenvalues[naive_ind]:.6f})"
            )
        follow_eigvec = eigenvectors.T[ind]
        step_nu = follow_eigvec.copy()
        nu = step_nu[-1]
        self.log(f"\tnu_{verbose}={nu:.8e}")
        # Scale eigenvector so that its last element equals 1. The
        # final is step is the scaled eigenvector without the last element.
        step = step_nu[:-1] / nu
        eigval = eigenvalues[ind]
        self.log(f"\teigenvalue_{verbose}={eigval:.8e}")
        return step, eigval, nu, follow_eigvec

    def filter_small_eigvals(self, eigvals, eigvecs, mask=False):
        small_inds = np.abs(eigvals) < self.small_eigval_thresh
        eigvals = eigvals[~small_inds]
        eigvecs = eigvecs[:, ~small_inds]
        small_num = sum(small_inds)
        self.log(
            f"Found {small_num} small eigenvalues in Hessian. Removed "
            "corresponding eigenvalues and eigenvectors."
        )
        assert small_num <= 6, (
            "Expected at most 6 small eigenvalues in cartesian hessian "
            f"but found {small_num}!"
        )
        if mask:
            return eigvals, eigvecs, small_inds
        else:
            return eigvals, eigvecs

    def log_negative_eigenvalues(self, eigvals, pre_str=""):
        neg_inds = eigvals < -self.small_eigval_thresh
        neg_eigval_str = np.array2string(eigvals[neg_inds], precision=6)
        self.log(f"{pre_str}Hessian has {neg_inds.sum()} negative eigenvalue(s).")
        self.log(f"\t{neg_eigval_str}")

    def housekeeping(self):
        """Calculate gradient and energy. Update trust radius and hessian
        if needed. Return energy, gradient and hessian for the current cycle."""
        gradient = self.geometry.gradient
        energy = self.geometry.energy
        self.forces.append(-gradient)
        self.energies.append(energy)
        self.log(f"    Energy: {energy: >12.6f} au")
        self.log(f"norm(grad): {np.linalg.norm(gradient): >12.6f} au / bohr (rad)")
        self.log(f" rms(grad): {np.sqrt(np.mean(gradient**2)): >12.6f} au / bohr (rad)")

        can_update = (
            # Allows gradient differences
            len(self.forces) > 1
            and (self.forces[-2].shape == gradient.shape)
            and len(self.coords) > 1
            # Coordinates may have been rebuilt. Take care of that.
            and (self.coords[-2].shape == self.coords[1].shape)
            and len(self.energies) > 1
        )
        if can_update:
            if self.trust_update:
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

        resetted = not can_update
        return energy, gradient, H, eigvals, eigvecs, resetted

    def get_augmented_hessian(self, eigvals, gradient, alpha=1.0):
        dim_ = eigvals.size + 1
        H_aug = np.zeros((dim_, dim_))
        H_aug[: dim_ - 1, : dim_ - 1] = np.diag(eigvals / alpha)
        H_aug[-1, :-1] = gradient
        H_aug[:-1, -1] = gradient

        H_aug[:-1, -1] /= alpha

        return H_aug

    def get_alpha_step(self, cur_alpha, rfo_eigval, step_norm, eigvals, gradient):
        # Derivative of the squared step w.r.t. alpha
        numer = gradient ** 2
        denom = (eigvals - rfo_eigval * cur_alpha) ** 3
        quot = np.sum(numer / denom)
        self.log(f"quot={quot:.6f}")
        dstep2_dalpha = 2 * rfo_eigval / (1 + step_norm ** 2 * cur_alpha) * quot
        self.log(f"analytic deriv.={dstep2_dalpha:.6f}")
        # Update alpha
        alpha_step = (
            2 * (self.trust_radius * step_norm - step_norm ** 2) / dstep2_dalpha
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
            rfo_step_, eigval_min, nu, self.prev_eigvec_min = self.solve_rfo(
                H_aug, "min", prev_eigvec=self.prev_eigvec_min
            )
            rfo_norm_ = np.linalg.norm(rfo_step_)
            self.log(f"norm(rfo step)={rfo_norm_:.6f}")

            if (rfo_norm_ < self.trust_radius) or abs(
                rfo_norm_ - self.trust_radius
            ) <= 1e-3:
                step_ = rfo_step_
                break

            alpha_step = self.get_alpha_step(
                alpha, eigval_min, rfo_norm_, eigvals, gradient_
            )
            alpha += alpha_step
            self.log("")
        # Otherwise, use trust region newton step
        else:
            self.log(
                "RS algorithm did not produce a desired step length in "
                f"{self.max_micro_cycles} micro cycles. Trying RFO with α=1.0."
            )
            H_aug = self.get_augmented_hessian(eigvals, gradient_, alpha=1.0)
            rfo_step_, eigval_min, nu, _ = self.solve_rfo(H_aug, "min")
            rfo_norm_ = np.linalg.norm(rfo_step_)

            # This should always be True if the above algorithm failed but we
            # keep this line nonetheless,  to make it more obvious.
            if rfo_norm_ > self.trust_radius:
                self.log(
                    f"Proposed RFO step with norm {rfo_norm_:.4f} is outside trust "
                    f"radius Δ={self.trust_radius:.4f}. "
                )
                step_ = self.get_newton_step_on_trust(
                    eigvals, eigvecs, gradient, transform=False
                )
                # Simple, downscaled RFO step
                # step_ = rfo_step_ / rfo_norm_ * self.trust_radius
            else:
                step_ = rfo_step_

        # Transform step back to original basis
        step = eigvecs.dot(step_)
        return step

    @staticmethod
    def get_shifted_step_trans(eigvals, gradient_trans, shift):
        return -gradient_trans / (eigvals + shift)

    @staticmethod
    def get_newton_step(eigvals, eigvecs, gradient):
        return eigvecs.dot(eigvecs.T.dot(gradient) / eigvals)

    def get_newton_step_on_trust(self, eigvals, eigvecs, gradient, transform=True):
        """Step on trust-radius.

        See Nocedal 4.3 Iterative solutions of the subproblem
        """
        min_ind = eigvals.argmin()
        min_eigval = eigvals[min_ind]
        pos_definite = (eigvals > 0.0).all()
        gradient_trans = eigvecs.T.dot(gradient)

        # This will be also be True when we come close to a minimizer,
        # but then the Hessian will also be positive definite and a
        # simple Newton step will be used.
        hard_case = abs(gradient_trans[min_ind]) <= 1e-6
        self.log(f"Smallest eigenvalue: {min_eigval:.6f}")
        self.log(f"Positive definite Hessian: {pos_definite}")
        self.log(f"Hard case: {hard_case}")

        def get_step(shift):
            return -gradient_trans / (eigvals + shift)

        # Unshifted Newton step
        newton_step_trans = get_step(0.0)
        newton_norm = np.linalg.norm(newton_step_trans)

        def on_trust_radius_lin(step):
            return 1 / self.trust_radius - 1 / np.linalg.norm(step)

        def finalize_step(shift):
            step = get_step(shift)
            if transform:
                step = eigvecs.dot(step)
            return step

        # Simplest case. Positive definite Hessian and predicted step is
        # already in trust radius.
        if pos_definite and newton_norm <= self.trust_radius:
            self.log("Using unshifted Newton step.")
            return eigvecs.dot(newton_step_trans)

        # If the Hessian is not positive definite or if the step is too
        # long we have to determine the shift parameter lambda.
        rs_kwargs = {
            "f": lambda shift: on_trust_radius_lin(get_step(shift)),
            "xtol": 1e-3,
            # Would otherwise be chosen automatically, but we set it
            # here explicitly for verbosity.
            "method": "brentq",
        }

        def root_search(bracket):
            rs_kwargs.update(
                {
                    "bracket": bracket,
                    "x0": bracket[0] + 1e-3,
                }
            )
            res = root_scalar(**rs_kwargs)
            return res

        BRACKET_END = 1e10
        if not hard_case:
            bracket_start = 0.0 if pos_definite else -min_eigval + 1e-2
            bracket = (bracket_start, BRACKET_END)
            try:
                res = root_search(bracket)
                assert res.converged
                return finalize_step(res.root)
            # ValueError may be raised when the function values for the
            # initial bracket have the same sign. If so, we continue with
            # treating it as a hard case.
            except ValueError:
                pass

        # Now we would try the bracket (-b2, -b1). The resulting step should have
        # a suitable length, but the (shifted) Hessian would have an incorrect
        # eigenvalue spectrum (not positive definite). To solve this we use a
        # different formula to calculate the step.
        mask = np.ones_like(gradient_trans)
        mask[min_ind] = 0
        mask = mask.astype(bool)
        without_min = gradient_trans[mask] / (eigvals[mask] - min_eigval)
        try:
            tau = sqrt(self.trust_radius ** 2 - (without_min ** 2).sum())
            step_trans = [tau] + (-without_min).tolist()
        # Hard case. Search in open interval (endpoints not included)
        # (-min_eigval, inf).
        except ValueError:
            bracket = (-min_eigval + 1e-6, BRACKET_END)
            res = root_search(bracket)
            if res.converged:
                return finalize_step(res.root)

        if not transform:
            return step_trans

        return eigvecs.dot(step_trans)

    @staticmethod
    def quadratic_model(gradient, hessian, step):
        return step.dot(gradient) + 0.5 * step.dot(hessian).dot(step)

    @staticmethod
    def rfo_model(gradient, hessian, step):
        return HessianOptimizer.quadratic_model(gradient, hessian, step) / (
            1 + step.dot(step)
        )

    def get_step_func(self, eigvals, gradient, grad_rms_thresh=1e-2):
        positive_definite = (eigvals < 0).sum() == 0
        gradient_small = rms(gradient) < grad_rms_thresh

        if self.adapt_step_func and gradient_small and positive_definite:
            return self.get_newton_step_on_trust, self.quadratic_model
        # RFO fallback
        else:
            return self.get_rs_step, self.rfo_model
