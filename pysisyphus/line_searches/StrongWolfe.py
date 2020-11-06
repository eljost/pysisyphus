# [1] Nocedal, Numerical Optimization

from pysisyphus.line_searches.LineSearch import LineSearch, \
                                                LineSearchNotConverged


from pysisyphus.line_searches.interpol import interpol_alpha_quad, interpol_alpha_cubic


class StrongWolfe(LineSearch):

    def __init__(self, *args, alpha_max=10., fac=2, **kwargs):
        """Wolfe line search.

        Uses only energy & gradient evaluations.

        See [1], Chapter 3, Line Search methods, Section 3.5 p. 60."""

        kwargs["cond"] = "strong_wolfe"
        super().__init__(*args, **kwargs)

        self.alpha_max = float(alpha_max)
        self.fac = fac

    def zoom(self, alpha_lo, alpha_hi, phi_lo,
             phi_alpha_=None, alpha_0_=None, max_cycles=10):
        phi0, dphi0 = self.get_phi_dphi("fg", 0)

        alphas = list()
        phi_alphas = list()
        if phi_alpha_:
            phi_alphas = [phi_alpha_, ]
        if alpha_0_:
            alphas = [alpha_0_, ]

        for j in range(max_cycles):
            # Interpoaltion of alpha between alpha_lo, alpha_hi
            #
            # Try cubic interpolation if at least two additional alphas and
            # corresponding phi_alpha values are available beside alpha = 0.
            if len(phi_alphas) > 1:
                alpha_prev = alphas[-1]
                phi_alpha_prev = phi_alphas[-1]
                alpha_j = interpol_alpha_cubic(phi0, dphi0,
                                               phi_alpha_, phi_alpha_prev,
                                               alpha_0_, alpha_prev
                )
            # Try quadratic interpolation if only one additional alpha and
            # corresponding phi_alpha value are available beside alpha = 0.
            elif len(phi_alphas) == 1:
                alpha_j = interpol_alpha_quad(phi0, dphi0, phi_alpha_, alpha_0_)
            # Fallback to simple bisection
            else:
                alpha_j = (alpha_lo + alpha_hi) / 2

            phi_j = self.get_phi_dphi("f", alpha_j)
            # Store the values so they can be reused for cubic interpolation
            alphas.append(alpha_j)
            phi_alphas.append(phi_j)

            # True if alpha is still too big or if the function value
            # increased compared to the previous cycle.
            if not self.sufficiently_decreased(alpha_j) or (phi_j > phi_lo):
                # Shrink interval to (alpha_lo, alpha_j)
                alpha_hi = alpha_j
                continue

            # If line search converged LineSearchConverged will be raised in
            # this call.
            dphi_j = self.get_phi_dphi("g", alpha_j)

            if (dphi_j * (alpha_hi - alpha_lo)) >= 0:
                alpha_hi = alpha_lo
            # Shrink interval to (alpha_j, alpha_hi)
            alpha_lo = alpha_j

    def run_line_search(self):
        phi0, dphi0 = self.get_phi_dphi("fg", 0)

        alpha_prev = 0
        phi_prev = phi0
        if self.alpha_init is not None:
            alpha_i = self.alpha_init
        else:
            alpha_i = 1.0

        for i in range(self.max_cycles):
            phi_i = self.get_phi_dphi("f", alpha_i)
            phi_rose = (phi_i >= phi_prev)
            # In [1] this condition is given with (if not sufficiently_decreased ...)
            # I guess this may give problems if the initial step is too small; then
            # suff. decr. is also not fullfilled ...
            if not self.sufficiently_decreased(alpha_i) or (phi_rose and i > 0):
                return self.zoom(alpha_prev, alpha_i, phi_prev, phi_i, alpha_i)

            dphi_i = self.get_phi_dphi("g", alpha_i)
            if self.strong_curvature_condition(alpha_i):
                break

            if dphi_i >= 0:
                return self.zoom(alpha_i, alpha_prev, phi_i, phi_alpha_=phi_i, alpha_0_=alpha_i)
            alpha_prev = alpha_i
            alpha_i = min(self.fac * alpha_i, self.alpha_max)
        # Premature abort of the loop may happen through LineSearchConverged being raised
        else:
            raise LineSearchNotConverged

        return alpha_i
