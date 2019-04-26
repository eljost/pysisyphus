#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
#     [4] 10.1007/s002140050387
#         Bofill, 1998, Restricted-Step-RFO
#     [5] https://link.springer.com/article/10.1007/s00214-016-1847-3
#         Birkholz, 2016

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.hessian_updates import (bfgs_update,
                                                   flowchart_update,
                                                   mod_flowchart_update)


class RSRFOptimizer(Optimizer):
    """Optimizer to find first-order saddle points."""

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }
    update_dict = {
        "flowchart": flowchart_update,
        "mod_flowchart": mod_flowchart_update,
        "bfgs": bfgs_update,
    }

    def __init__(self, geometry, trust_radius=.1, calc_hess=None,
                 recalc_hess=None,
                 max_micro_cycles=50, hess_update="flowchart",
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius0 = trust_radius
        self.calc_hess = calc_hess
        self.recalc_hess = recalc_hess
        self.max_micro_cycles = max_micro_cycles
        self.hess_update = hess_update
        self.update_func = self.update_dict[self.hess_update]
        # Trust radius thresholds
        self.min_trust_radius = 0.25 * self.trust_radius0
        # self.max_trust_radius = 2 * self.trust_radius0
        self.max_trust_radius = 2

        self.alpha0 = 1
        self.alpha_max = 1e8
        self.trust_radius = self.trust_radius0
        self.H = None
        self.predicted_energy_changes = list()

    def prepare_opt(self):
        if self.calc_hess:
            self.H = self.geometry.hessian
        else:
            self.H = self.geometry.get_initial_hessian()

    def update_trust_radius(self, coeff, last_step_norm):
        # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1
        # The trust radius update proposed in [4], Sec. 3.1.5 seems
        # a bit crazy.
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius/4,
                                    self.min_trust_radius)
            self.log("Decreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
            self.trust_radius = min(self.trust_radius*2,
                                    self.max_trust_radius)
            self.log("Increasing trust radius.")
        else:
            self.log("Keeping current trust radius")
            return
        self.log(f"Updated trust radius: {self.trust_radius:.6f}")

    def solve_rfo(self, rfo_mat, kind="min"):
        # So if I use eig instead of eigh here it even works ...
        # my bad, ahhh! The unscaled RFO matrix may be symmetric,
        # but the scaled ones aren't anymore.
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

    def optimize(self):
        forces = self.geometry.forces
        fn = np.linalg.norm(forces)
        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        # Hessian recalculation/update
        if (self.recalc_hess and (self.cur_cycle > 0)
            and (self.cur_cycle % self.recalc_hess) == 0):
            self.log("Recalculating exact hessian")
            self.H = self.geometry.hessian
        elif self.cur_cycle > 0:
            # Gradient difference
            dg = -(self.forces[-1] - self.forces[-2])
            # Coordinate difference
            dx = self.steps[-1]
            H_update, key = self.update_func(self.H, dx, dg)
            self.H += H_update
            self.log(f"{key} hessian update")

        if self.cur_cycle > 0:
            actual_energy_change = self.energies[-1] - self.energies[-2]
            self.log(f"actual energy change: {actual_energy_change:.4e}")
            predicted_energy_change = self.predicted_energy_changes[-1]
            self.log(f"predicted energy change: {predicted_energy_change:.4e}")
            coeff = actual_energy_change / predicted_energy_change
            self.log(f"energy change, actual/predicted={coeff:.6f}")
            last_step_norm = np.linalg.norm(self.steps[-1])
            self.log(f"Current trust radius is {self.trust_radius}")
            self.update_trust_radius(coeff, last_step_norm)

        H = self.H
        if self.geometry.internal:
            H = self.geometry.internal.project_hessian(self.H)
        eigvals, eigvecs = np.linalg.eigh(H)

        # Transform to eigensystem of hessian
        forces_trans = eigvecs.T.dot(forces)

        # Minimize energy along all modes
        min_mat = np.asarray(np.bmat((
            (np.diag(eigvals), -forces_trans[:,None]),
            (-forces_trans[None,:], [[0]])
        )))

        alpha = self.alpha0
        min_diag_indices = np.diag_indices(eigvals.size)
        for mu in range(self.max_micro_cycles):
            assert alpha > 0, "alpha should not be negative"
            self.log(f"RS-RFO micro cycle {mu:02d}, alpha={alpha:.6f}")
            # We only have to update one eigenvalue
            min_mat_scaled = min_mat.copy()
            min_mat_scaled[min_diag_indices] /= alpha
            min_mat_scaled[:-1,-1] /= alpha
            rfo_step, eigval_min, nu_min = self.solve_rfo(min_mat_scaled, "min")

            # As of Eq. (8a) of [4] max_eigval and min_eigval also
            # correspond to:
            # eigval_min_ = -forces_trans.dot(rfo_step)
            # np.testing.assert_allclose(eigval_min, eigval_min_)

            # Create the full PRFO step
            rfo_norm = np.linalg.norm(rfo_step)
            self.log(f"rfo_norm={rfo_norm:.6f}")

            inside_trust = rfo_norm < self.trust_radius + 1e-3
            if inside_trust:
                self.log("step is inside trust radius. breaking.")
                break
            elif alpha > self.alpha_max:
                print("alpha > alpha_max. breaking.")
                break

            # Derivative of the squared step w.r.t. alpha
            tval = 2*eigval_min/(1+rfo_norm**2 * alpha)
            numer = forces_trans**2
            denom = (eigvals - eigval_min * alpha)**3
            quot = np.sum(numer / denom)
            self.log(f"quot={quot:.6f}")
            dstep2_dalpha = (2*eigval_min/(1+rfo_norm**2 * alpha)
                             * np.sum(forces_trans**2
                                      / ((eigvals - eigval_min * alpha)**3)
                               )
            )
            self.log(f"analytic deriv.={dstep2_dalpha:.6f}")
            # Update alpha
            alpha_step = (2*(self.trust_radius*rfo_norm - rfo_norm**2)
                          / dstep2_dalpha
            )
            self.log(f"alpha_step={alpha_step:.4f}")
            alpha += alpha_step
            self.log("")

        predicted_energy_change = 1/2 * eigval_min / nu_min**2
        self.predicted_energy_changes.append(predicted_energy_change)
        self.log(f"predicted_energy_change={predicted_energy_change:.6e}")
        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(rfo_step)
        step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={step_norm:.6f}")

        self.log("")
        return step
