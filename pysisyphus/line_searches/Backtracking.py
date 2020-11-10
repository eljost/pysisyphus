from pysisyphus.line_searches.LineSearch import (
    LineSearch,
    LineSearchNotConverged,
)

from pysisyphus.line_searches.interpol import interpol_alpha_quad, interpol_alpha_cubic
from pysisyphus.optimizers.poly_fit import cubic_fit, quartic_fit


class Backtracking(LineSearch):
    def __init__(self, *args, rho_lo=5e-2, rho_hi=0.9, use_grad=False, **kwargs):
        """Backtracking line search enforcing Armijo conditions.

        Uses only energy evaluations.

        See [1], Chapter 3, Line Search methods, Section 3.1 p. 31 and
        Section 3.5 p. 56."""

        kwargs["cond"] = "armijo"
        super().__init__(*args, **kwargs)

        self.rho_lo = float(rho_lo)
        self.rho_hi = float(rho_hi)
        self.use_grad = use_grad

    def alpha_new_from_phi(self, cycle, phi0, dphi0, alpha, alpha_prev):
        phi_i = self.get_phi_dphi("f", alpha)
        self.log(f"\tCycle {cycle:02d}: alpha={alpha:.6f}, ϕ={phi_i:.6f} au")

        if cycle == 0:
            # Quadratic interpolation
            alpha_new = interpol_alpha_quad(phi0, dphi0, phi_i, alpha)
            type_ = "quadratic"
        else:
            # Cubic interpolation
            phi_prev = self.get_phi_dphi("f", alpha_prev)
            alpha_new = interpol_alpha_cubic(
                phi0, dphi0, phi_prev, phi_i, alpha_prev, alpha
            )
            type_ = "cubic"
        return alpha_new, type_

    def alpha_new_from_phi_dphi(self, cycle, phi0, dphi0, alpha):
        phi_i, dphi_i = self.get_phi_dphi("fg", alpha)
        self.log(f"\tCycle {cycle:02d}: α={alpha:.6f}, ϕ={phi_i:.6f} au")

        # First we try a constrained quartic polynomial
        res = quartic_fit(phi0, phi_i, dphi0, dphi_i)
        type_ = "quartic"
        # If the quartic poly failed, we continue with a cubic polynomial
        if res is None:
            res = cubic_fit(phi0, phi_i, dphi0, dphi_i)
            type_ = "cubic"
        # If the cubic poly failed we resort to bisection. ohoh
        if res is None:
            alpha_new = 0.5 * alpha
            type_ = "bisection"
        else:
            alpha_new = res.x * alpha
        return alpha_new, type_

    def run_line_search(self):
        phi0, dphi0 = self.get_phi_dphi("fg", 0)

        alpha_prev = None
        alpha = self.alpha_init
        for i in range(self.max_cycles):
            if self.use_grad:
                alpha_new, type_ = self.alpha_new_from_phi_dphi(i, phi0, dphi0, alpha)
            else:
                alpha_new, type_ = self.alpha_new_from_phi(i, phi0, dphi0, alpha, alpha_prev)
            self.log(f"\tNew α from {type_} interpolation: {alpha_new:.6f}")

            lower_bound = alpha * self.rho_lo
            upper_bound = alpha * self.rho_hi
            if alpha_new < lower_bound:
                self.log("\tNew α is too small!")
            if alpha_new > upper_bound:
                self.log("\tNew α is too big!")

            # Assert that alpha doesn't change too much compared to the previous alpha
            alpha_new = min(alpha_new, upper_bound)
            alpha_new = max(alpha_new, lower_bound)
            alpha_prev = alpha
            alpha = alpha_new
            self.log(f"\tNext α: {alpha:.6f}\n")

        raise LineSearchNotConverged
