#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np
from scipy.optimize import minimize

from pysisyphus.optimizers.Optimizer import Optimizer


class RFOptimizer(Optimizer):

    def __init__(self, geometry, trust_radius=0.3, gdiis_thresh=-1, **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius
        self.gdiis_thresh = gdiis_thresh

        self.min_trust_radius = 0.25*self.trust_radius
        self.trust_radius_max = 5*self.trust_radius
        self.predicted_energy_changes = list()
        self.actual_energy_changes = list()

        # self.hessians = list()
        self.trust_radii = list()
        self.rfo_steps = list()

        self.error_vecs = list()
        # Allow greaeter deviations from the reference step when
        # more error vectors are used.
        self.gdiis_cos_cutoffs = {
            2: 0.97,
            3: 0.84,
            2: 0.80,
            3: 0.75,
            4: 0.71,
            5: 0.67,
            6: 0.62,
            7: 0.56,
            8: 0.49,
            9: 0.41,
        }

    def prepare_opt(self):
        self.H = self.geometry.get_initial_hessian()

    def keep(self):
        # self.hessians.append(self.H.copy())
        self.trust_radii.append(self.trust_radius)
        # self.log("!Saving hessian every iteration!")

    def bfgs_update(self):
        # Eq. (44) in [1]
        dx = self.coords[-1] - self.coords[-2]
        dg = -(self.forces[-1] - self.forces[-2])
        second_term = np.outer(dg, dg) / np.inner(dg, dx)
        third_term = (self.H.dot(np.outer(dx, dx)).dot(self.H)
                      / dx.dot(self.H.dot(dx)))
        self.H += second_term - third_term

    def quadratic_approx(self, step):
        E0 = self.energies[-1]
        g = -self.forces[-1]
        return E0 + np.inner(g, step) + 0.5*step.dot(self.H.dot(step))

    def predicted_change(self, step):
        E0 = self.energies[-1]
        g = -self.forces[-1]
        return np.inner(g, step) + 0.5*step.dot(self.H.dot(step))

    def find_step(self, step_guess):
        ineq_fun = lambda step: self.trust_radius - np.linalg.norm(step)
        constr = {
                "type": "ineq",
                "fun": ineq_fun,
        }
        res = minimize(self.quadratic_approx,
                 step_guess,
                 method="SLSQP",
                 constraints=constr)
        if not res.success:
            self.log("LQA optimization failed!")
        step = res.x
        self.log(f"LQA minimum: {res.fun:.6f} au")
        self.log(f"Optimized step norm: {np.linalg.norm(step):.4f}")
        original_direction = step_guess / np.linalg.norm(step_guess)
        new_direction = step / np.linalg.norm(step)
        dir_rad = np.arccos(original_direction.dot(new_direction))
        dir_deg = np.rad2deg(dir_rad)
        self.log(f"Angle between original and new direction: {dir_deg:.2f}°")
        return step

    def update_trust_radius(self):
        # [3] Chapter 4, Algorithm 4.1
        actual_reduction = self.energies[-2] - self.energies[-1]
        predicted = self.predicted_energy_changes[-1]
        actual = self.actual_energy_changes[-1]
        reduction_ratio = actual / predicted
        self.log(f"Predicted energy reduction: {predicted:.4e} au")
        self.log(f"Actual energy reduction: {actual:.4e} au")
        self.log(f"Reduction ratio: {reduction_ratio:.5f}")
        last_step_norm = np.linalg.norm(self.steps[-1])
        if reduction_ratio < 0.25:
            self.trust_radius = max(0.25*self.trust_radius,
                                    self.min_trust_radius)
            self.log("Decreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [4], Appendix, step size and direction control
        elif reduction_ratio > 0.5 and (last_step_norm >= .8*self.trust_radius):
            self.trust_radius = min(2*self.trust_radius, self.trust_radius_max)
            self.log("Increasing trust radius.")
        self.log(f"Trust radius: {self.trust_radius:.3f}")

    def gdiis(self, ref_step, max_err_vecs=5):
        np.set_printoptions(suppress=True, precision=4)
        # max of 5 diis vectors as proposed by swart 2006
        max_vecs = min(max_err_vecs, len(self.error_vecs))
        valid_diis_step = None
        valid_A = None
        ref_step_dir = ref_step / np.linalg.norm(ref_step)
        for i in range(max_vecs):#len(self.error_vecs)):
            inds = list(range(i+1))
            if i == 0:
                continue
            A = np.zeros((i+1, i+1))
            # Start with the latest point and add previous points one by one
            err_vecs = np.array(self.error_vecs[::-1])[inds]
            # Scale error vector so that norm of smallest error vector is 1
            err_norms = np.linalg.norm(err_vecs, axis=1)
            # print("unscaled error vector norms", err_norms)
            scale_factor = 1 / err_norms.min()
            err_vecs  *= scale_factor
            err_norms = np.linalg.norm(err_vecs, axis=1)
            # print("scaled error vector norms", err_norms)
            for i, e1 in enumerate(err_vecs):
                for j in range(i, len(err_vecs)):
                    e2 = err_vecs[j]
                    A[i, j] = e1.dot(e2)
                    A[j, i] = A[i, j]

            vec_num = len(err_vecs)
            A_ = np.ones((vec_num+1, vec_num+1))
            for i, g1 in enumerate(err_vecs):
                for j, g2 in enumerate(err_vecs):
                    A_[i, j] = g1.dot(g2)
            A_[-1, -1] = 0
            b = np.zeros(vec_num+1)
            b[-1] = 1
            cs_lambda = np.linalg.solve(A_, b)
            *cs_, lam_ = cs_lambda


            # print(f"gdiis, using {len(err_vecs)} error vector(s).")
            # print("A")
            # print(A)
            cr = np.linalg.solve(A, np.ones(A.shape[0]))
            # Rule 1: Check for c/|r²| elements > 1e8
            # print("cr", cr, "sum", np.sum(cr))
            if any(np.abs(cr) > 1e8):
                # print("Found big cr element!")
                break
            cs = cr / np.sum(cr)
            np.testing.assert_allclose(cs_, cs)
            # print("normalized c", cs)
            # Rule 2: Absolute value of sum of pos. (neg.) coefficients < 15
            pos_c_sum = np.sum(cs[cs > 0])
            neg_c_sum = np.sum(cs[cs < 0])
            if (pos_c_sum > 15) or (np.abs(neg_c_sum) > 15):
                # print("Sum of coefficients > 15!")
                break
            # print("sum pos c", pos_c_sum)
            # print("sum neg c", neg_c_sum)
            last_coords = np.array(self.coords[::-1])[inds]
            last_steps = np.array([ref_step] + self.steps[::-1])[inds]
            diis_coords = np.sum(cs[:,None]*last_coords, axis=0)
            diis_step = np.sum(cs[:,None]*last_steps, axis=0)
            # _ = self.coords[-1] - diis_coords
            # diis_step = self.coords[-1] - diis_coords
            # import pdb; pdb.set_trace()
            # Rule 3: Compare magnitude of reference and DIIS step and
            # allow only 10 times the magnitude of the reference step
            # for the DIIS step.
            diis_step_norm = np.linalg.norm(diis_step)
            ref_step_norm = np.linalg.norm(ref_step)
            if diis_step_norm > 10*ref_step_norm:
                factor = (10*ref_step_norm) / diis_step_norm
                diis_step *= factor
                print("scaled down diis step")

            # Rule 4: Similar direction as reference step
            diis_step_dir =  diis_step / np.linalg.norm(diis_step)
            cos_cutoff = self.gdiis_cos_cutoffs[cs.size]
            cos_ = diis_step_dir.dot(ref_step_dir)
            if cos_ < cos_cutoff:
                print("cos_=", cos_, "cutoff", cos_cutoff)
                # print("directions of reference and diis differ too much. breaking!")
                break
            valid_diis_step = diis_step
            valid_A = A
        if valid_diis_step is not None:
            print(f"Did GDIIS with {valid_A.shape[0]} error vectors.")
        return valid_diis_step

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            # Predicted changes
            # predicted_energy_change = self.quadratic_approx(self.steps[-1])
            predicted_energy_change = self.predicted_change(self.steps[-1])
            self.predicted_energy_changes.append(predicted_energy_change)
            # Actual change
            self.log(f"Last two energies: {self.energies[-2:]}")
            actual_energy_change = self.energies[-1] - self.energies[-2]
            self.actual_energy_changes.append(actual_energy_change)
            self.update_trust_radius()
            self.bfgs_update()

        if self.geometry.internal:
            self.H = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it.
            self.H = (self.H + self.H.T) / 2

        # Eq. (56) in [1]
        aug_hess = np.bmat(
                    ((self.H, gradient[:,None]),
                     (gradient[None,:], [[0]]))
        )
        eigvals, eigvecs = np.linalg.eigh(aug_hess)
        self.log(f"First 5 eigenvalues: {np.array2string(eigvals[:5], precision=2, suppress_small=True)}")
        self.log(f"Number of negative eigenvalues: {eigvals[eigvals < 0].size}")
        self.log(f"Lowest eigenvalue: {eigvals[0]:.6f}")
        # Select eigenvector corresponding to smallest eigenvalue.
        # As the eigenvalues are sorted in ascending order eigvals.argmin()
        # should always give 0...
        assert eigvals.argmin() == 0
        aug_step = eigvecs[:,0]
        # aug_step is currently a matrix. Convert it to an array.
        aug_step = np.asarray(aug_step).flatten()
        # Scale aug_step so the last element equals 1
        lambda_ = aug_step[-1]
        self.log(f"lambda: {lambda_:.6f}")
        aug_step /= lambda_
        step = aug_step[:-1]

        self.keep()
        self.rfo_steps.append(step)

        # Restrict elements of the the step vector to an allowed maximum
        # if they exceed it.
        # step[np.abs(step) > 0.3] = 0.3
        step_norm = np.linalg.norm(step)
        self.log(f"Unscaled norm(step): {step_norm:.4f}")
        # We use a trust region method instead
        if step_norm > self.trust_radius:
            step = self.find_step(step)
        self.log(f"norm(step): {np.linalg.norm(step):.4f}")
        self.log("")

        self.error_vecs.append(step)
        rms_rfo_step = np.sqrt(np.mean(step**2))
        # Start GDIIS
        if rms_rfo_step < self.gdiis_thresh:
            diis_step = self.gdiis(step)
            # print(diis_step)
            if diis_step is not None:
                print("found valid diis step")
                step = diis_step

        return step
