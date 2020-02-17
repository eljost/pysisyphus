import numpy as np

from pysisyphus.line_searches.LineSearch import LineSearch, \
                                                LineSearchNotConverged


class HagerZhang(LineSearch):

    def __init__(self, *args, alpha_prev=None, f_prev=None, dphi0_prev=None,
                 quad_step=False, eps=1e-6, theta=0.5, gamma=0.5, rho=5,
                 psi_0=.01, psi_1=.1, psi_2=2., psi_low=0.1, psi_hi=10,
                 Delta=.7, omega=1e-3, max_bisects=10, **kwargs):

        kwargs["cond"] = "wolfe"
        super().__init__(*args, **kwargs)

        self.alpha_prev = alpha_prev
        self.f_prev = f_prev
        self.dphi0_prev = dphi0_prev
        self.quad_step = quad_step

        self.eps = eps
        self.theta = theta
        self.gamma = gamma
        self.rho = rho
        self.psi_0 = psi_0
        self.psi_1 = psi_1
        self.psi_2 = psi_2
        self.psi_low = psi_low
        self.psi_hi = psi_hi
        self.Delta = Delta
        self.omega = omega
        self.max_bisects = max_bisects

    def prepare_line_search(self):
        super().prepare_line_search()

        self.epsk = self.eps * abs(self.get_fg("f", 0.))

    def bisect(self, a, b):
        """Bisect interval [a, b]."""
        for i in range(self.max_bisects):
            # U3 a.
            d = (1 - self.theta)*a + self.theta*b
            dphi_d = self.get_phi_dphi("g", d)
            if dphi_d >= 0:
                return a, d

            phi_d = self.get_phi_dphi("f", d)
            # U3 b.
            # If (dphi_d > 0) we would already have returned above...
            if phi_d <= self.phi0 + self.epsk:
                a = d
            # U3 c.
            elif phi_d > self.phi0 + self.epsk:
                b = d
        raise Exception("Bisect failed!")

    def interval_update(self, a, b, c):
        """Narrows down the bracketing interval."""
        # U0
        if not (a < c < b):
            return a, b

        phi_c, dphi_c = self.get_phi_dphi("fg", c)
        # U1, sign of slope projection changed. We already passed the minimum.
        if dphi_c >= 0:
            return a, c
        # U2, we are moving towards the minimum.
        elif phi_c <= self.phi0 + self.epsk:
            return c, b

        # U3, phi_c increased above phi0, so we probably passed the minimum.
        return self.bisect(a, c)

    def secant(self, a, b):
        """Take secant step."""
        dphia = self.get_phi_dphi("g", a)
        dphib = self.get_phi_dphi("g", b)
        return (a*dphib - b*dphia) / (dphib - dphia)

    def double_secant(self, a, b):
        """Take secant² step."""
        c = self.secant(a, b)
        A, B = self.interval_update(a, b, c)
        cB_close = np.isclose(c, B)
        cA_close = np.isclose(c, A)

        if cB_close:
            c_dash = self.secant(b, B)
        elif cA_close:
            c_dash = self.secant(a, A)

        if cB_close or cA_close:
            a_dash, b_dash = self.interval_update(A, B, c_dash)
        else:
            a_dash, b_dash = A, B
        return a_dash, b_dash

    def bracket(self, c):
        """Generate initial interval [a, b] that satisfies the opposite
        slope condition (dphi(a) < 0, dphi(b) > 0).
        """
        cs = list()
        p0epsk = self.phi0 + self.epsk
        for j in range(10):
            cs.append(c)

            dphi_j = self.get_phi_dphi("g", c)

            if (dphi_j >= 0) and (j == 0):
                return 0, c

            phi_j = self.get_phi_dphi("f", c)
            if dphi_j >= 0:
                phi_inds = np.array([self.get_fg("f", c) for c in cs[:-1]]) <= p0epsk
                # See https://stackoverflow.com/a/8768734
                ci = len(phi_inds) - phi_inds[::-1].argmax() - 1
                return cs[ci], c
            elif phi_j > p0epsk:
                return self.bisect(0, c)

            c *= self.rho

    def norm_inf(self, arr):
        """Returns infinity norm of given array."""
        return np.linalg.norm(arr, np.inf)

    def initial(self):
        """Get an initial guess for alpha."""
        if (~np.isclose(self.x0, np.zeros_like(self.x0))).any():
            c = self.psi_0 * self.norm_inf(self.x0)/self.norm_inf(self.g0)
        elif not np.isclose(self.f0, 0):
            c = self.psi_0 * self.f0 / self.norm_inf(self.g0)**2
        else:
            c = 1
        return c

    def take_quad_step(self, alpha, g0_):
        """Try to get alpha for minimum step from quadratic interpolation."""
        fact = max(self.psi_low, g0_/(self.dphi0*self.psi_2))
        alpha_ = min(fact, self.psi_hi) * alpha
        phi_ = self.get_phi_dphi("f", alpha_)
        denom = 2*((phi_-self.phi0)/alpha_ - self.dphi0)
        f_temp = self.get_fg("f", alpha_)
        if denom > 0.:
            c = -self.dphi0*alpha_ / denom
            if f_temp > self.get_fg("f", 0):
                c = max(c, alpha_*1e-10)
        else:
            c = alpha
        return c

    def run_line_search(self):
        if self.alpha_init is None and self.alpha_prev:
            alpha_init = self.alpha_prev
        elif self.alpha_init is None and self.alpha_prev is None:
            alpha_init = self.initial()
        else:
            alpha_init = self.alpha_init

        if self.quad_step:
            g0_ = -2*abs(self.get_fg("f", 0)/alpha_init) if (self.dphi0_prev is None) \
                  else self.dphi0_prev
            alpha_init = self.take_quad_step(self.psi_2*alpha_init, g0_)
        # This may raise LineSearchConverged
        _ = self.get_phi_dphi("fg", alpha_init)

        # TODO: cubic interpolation for better alpha_init
        ak, bk = self.bracket(alpha_init)
        for k in range(self.max_cycles):
            if self.cond_func(ak):
                break
            # secant² step
            a, b = self.double_secant(ak, bk)
            if (b - a) > self.gamma*(bk - ak):
                # Bisection step
                c = (a + b)/2
                a, b = self.interval_update(a, b, c)
            ak, bk = a, b
        else:
            raise LineSearchNotConverged

        return ak
