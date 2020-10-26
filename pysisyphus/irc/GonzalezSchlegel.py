import numpy as np
from scipy.optimize import newton

from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.TableFormatter import TableFormatter


# [1] An improved algorithm for reaction path following
# http://aip.scitation.org/doi/pdf/10.1063/1.456010
# [2] Extension to internal coordinates (not implemented)
# https://pubs.acs.org/doi/pdf/10.1021/j100377a021


class GonzalezSchlegel(IRC):

    def __init__(self, geometry, max_micro_cycles=20, micro_step_thresh=1e-3,
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.max_micro_cycles = max_micro_cycles
        self.micro_step_thresh = micro_step_thresh

        self.pivot_coords = list()
        self.micro_coords = list()
        self.eye = np.eye(self.geometry.coords.size)

        micro_header = "# |dx| |tangent|".split()
        micro_fmts = ["d", ".2E", ".3E"]
        self.micro_formatter = TableFormatter(micro_header, micro_fmts, 10)

    def micro_step(self):
        """Constrained optimization on a hypersphere."""

        gradient = self.mw_gradient
        gradient_diff = gradient - self.prev_grad
        coords_diff = self.mw_coords - self.prev_coords
        self.prev_grad = gradient
        # Without copy we would only store the reference...
        self.prev_coords = self.mw_coords.copy()

        dH, _ = bfgs_update(self.mw_hessian, coords_diff, gradient_diff)
        self.mw_hessian += dH
        eigvals, eigvecs = np.linalg.eigh(self.mw_hessian)

        constraint = (0.5 * self.step_length)**2
        big = np.abs(eigvals) > 1e-8
        big_eigvals = eigvals[big]
        big_eigvecs = eigvecs[:, big]
        grad_star = big_eigvecs.T.dot(gradient)
        displ_star = big_eigvecs.T.dot(self.displacement)

        def get_dx(lambda_):
            return -(grad_star - lambda_*displ_star) / (big_eigvals - lambda_)

        def on_sphere(lambda_):
            p = displ_star + get_dx(lambda_)
            return p.dot(p) - constraint

        # Initial guess for λ.
        # λ must be smaller then the smallest eigenvector
        lambda_0 = big_eigvals[0]
        lambda_0 *= 1.5 if (lambda_0 < 0) else 0.5
        # Find the root with scipy
        lambda_ = newton(on_sphere, lambda_0, maxiter=500)

        # Calculate dx from optimized lambda in basis of Hessian eigenvectors and
        # transform back to mass-weighted Cartesians.
        dx = big_eigvecs.dot(get_dx(lambda_))
        self.displacement += dx
        self.mw_coords += dx

        displ_norm = np.linalg.norm(self.displacement)
        tangent = gradient - gradient.dot(self.displacement)/displ_norm * gradient

        return dx, tangent

    def step(self):
        grad0 = self.mw_gradient
        grad0_norm = np.linalg.norm(grad0)
        # For the BFGS update in the first micro step we use the original
        # point and the initial guess to calculate gradient and
        # coordinate differences.
        self.prev_grad = grad0
        self.prev_coords = self.mw_coords

        # Take a step against the gradient to the pivot point x*_k+1.
        pivot_step = 0.5*self.step_length * -grad0/grad0_norm
        pivot_coords = self.mw_coords + pivot_step
        self.pivot_coords.append(pivot_coords)

        # Initial guess for x'_k+1 (full step from prev_coords, or another
        # half step from the pivot point)
        self.mw_coords = pivot_coords + pivot_step

        # Initial displacement p' from the pivot point
        self.displacement = pivot_step

        micro_coords_ = list()
        self.table.print(self.micro_formatter.header)
        for i in range(self.max_micro_cycles):
            try:
                dx, tangent = self.micro_step()
            except RuntimeError:
                print("Constrained search did not converge!")
                self.converged = True
                return
            micro_coords_.append(self.mw_coords)
            norm_dx = np.linalg.norm(dx)
            norm_tangent = np.linalg.norm(tangent)
            self.table.print(self.micro_formatter.line(i+1, norm_dx, norm_tangent))

            if (np.linalg.norm(dx) <= self.micro_step_thresh):
                break
        else:
            self.logger.warning("Max micro cycles exceeded!")

        # calc = self.geometry.calculator
        # pc = self.pivot_coords[0]
        # c = self.coords
        # xs = (c[0], pc[0])
        # ys = (c[1], pc[1])
        # calc.plot_coords(xs, ys, show=True)

        self.micro_coords.append(np.array(micro_coords_))

    def postprocess(self):
        self.pivot_coords = np.array(self.pivot_coords)
        self.micro_coords = np.array(self.micro_coords)