from pysisyphus.line_searches.LineSearch import LineSearch, \
                                                LineSearchNotConverged

from pysisyphus.line_searches.interpol import interpol_alpha_quad, interpol_alpha_cubic


class Backtracking(LineSearch):

    def __init__(self, *args, rho_lo=5e-2, rho_hi=0.9, **kwargs):
        """Backtracking line search enforcing Armijo conditions.

        Uses only energy evaluations.

        See [1], Chapter 3, Line Search methods, Section 3.1 p. 31 and
        Section 3.5 p. 56."""

        kwargs["cond"] = "armijo"
        super().__init__(*args, **kwargs)

        self.rho_lo = float(rho_lo)
        self.rho_hi = float(rho_hi)

    def run_line_search(self):
        self.log("Starting backtracking line search")
        phi0, dphi0 = self.get_phi_dphi("fg", 0)

        alpha_prev = None
        alpha = self.alpha_init
        for i in range(self.max_cycles):
            phi_i = self.get_phi_dphi("f", alpha)
            self.log(f"\tCycle {i:02d}: alpha={alpha:.6f}, ϕ={phi_i:.6f}")

            if self.cond_func(alpha):
                self.log(f"\tLine search converged after {i} cycles.")
                break

            if i == 0:
                # Quadratic interpolation
                alpha_new = interpol_alpha_quad(phi0, dphi0, phi_i, alpha)
                type_ = "Quadratic"
            else:
                # Cubic interpolation
                phi_prev = self.get_phi_dphi("f", alpha_prev)
                alpha_new = interpol_alpha_cubic(phi0, dphi0, phi_prev, phi_i, alpha_prev, alpha)
                type_ = "Cubic"
            self.log(f"\tNew alpha from {type_}: {alpha_new:.6f}")

            lower_bound = alpha * self.rho_lo
            upper_bound = alpha * self.rho_hi
            if alpha_new < lower_bound:
                self.log("\tNew alpha is too big!")
            if alpha_new > upper_bound:
                self.log("\tNew alpha is too high!")

            # Assert that alpha doesn't change too much compared to the previous alpha   
            alpha_new = min(alpha_new, upper_bound)
            alpha_new = max(alpha_new, lower_bound)
            alpha_prev = alpha
            alpha = alpha_new
            self.log(f"\tAlpha for next cycles: {alpha:.6f}\n")
        else:
            raise LineSearchNotConverged

        return alpha
