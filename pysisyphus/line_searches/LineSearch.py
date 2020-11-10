import abc
from collections import namedtuple
import logging

import numpy as np


class LineSearchConverged(Exception):
    def __init__(self, alpha):
        self.alpha = alpha


class LineSearchNotConverged(Exception):
    pass


LineSearchResult = namedtuple(
    "LineSearchResult",
    "converged alpha f_new g_new f_evals df_evals dphi0",
    # defaults=( None, None, None, None,   None,    None),
)


class LineSearch(metaclass=abc.ABCMeta):
    def __init__(
        self,
        p,
        cond="armijo",
        x0=None,
        geometry=None,
        f=None,
        df=None,
        alpha_init=None,
        f0=None,
        g0=None,
        c1=0.1,
        c2=0.9,
        max_cycles=10,
        alpha_min=1e-6,
    ):
        self.p = p
        self.geometry = geometry
        self.f = f
        self.df = df

        geometry_supplied = self.geometry is not None
        assert geometry_supplied or (
            x0 is not None
        ), "Supply either 'geometry' or the starting coordinates 'x0'!"
        assert geometry_supplied or (self.f and self.df), (
            "Supply either 'geometry' with a calculator or the two functions "
            "'f' and 'df' to calculate the energy and its gradient!"
        )

        self.x0 = x0
        if self.geometry:
            self.f = lambda coords: self.geometry.get_energy_at(coords)
            self.df = lambda coords: -self.geometry.get_energy_and_forces_at(coords)[
                "forces"
            ]
            self.x0 = self.geometry.coords.copy()

        self.alpha_init = alpha_init
        self.f0 = f0
        self.g0 = g0
        self.c1 = c1
        self.c2 = c2
        self.max_cycles = max_cycles
        self.alpha_min = alpha_min

        # Store calculated energies & gradients
        self.alpha_fs = {}
        self.alpha_dfs = {}
        self.f_evals = 0
        self.df_evals = 0

        self.dphis = {}

        self.cond_funcs = {
            "armijo": self.sufficiently_decreased,
            "wolfe": self.wolfe_condition,
            "strong_wolfe": self.strong_wolfe_condition,
        }
        self.can_eval_cond_funcs = {
            "armijo": lambda alpha: alpha in self.alpha_fs,
            "wolfe": self.got_alpha_phi_dphi,
            "strong_wolfe": self.got_alpha_phi_dphi,
        }
        self.cond_func = self.cond_funcs[cond]
        self.can_eval_cond_func = self.can_eval_cond_funcs[cond]

        self.logger = logging.getLogger("optimizer")

    def log(self, message):
        self.logger.debug(message)

    def prepare_line_search(self):
        if self.f0 is None:
            self.phi0 = self.get_phi_dphi("f", 0)
        else:
            self.phi0 = self.f0
            self.alpha_fs[0.0] = self.f0
        if self.g0 is None:
            self.dphi0 = self.get_phi_dphi("g", 0)
        else:
            self.dphi0 = self.g0.dot(self.p)
            self.alpha_dfs[0.0] = self.g0

    def check_alpha(self, alpha):
        if (alpha != 0.0) and (np.isnan(alpha) or (alpha < self.alpha_min)):
            raise LineSearchNotConverged()

    def _phi(self, alpha):
        alpha = float(alpha)
        self.check_alpha(alpha)
        try:
            f_alpha = self.alpha_fs[alpha]
        except KeyError:
            self.log(f"\tEvaluating energy for alpha={alpha:.6f}")
            f_alpha = self.f(self.x0 + alpha * self.p)
            self.f_evals += 1
            self.alpha_fs[alpha] = f_alpha
        return f_alpha

    def _dphi(self, alpha):
        alpha = float(alpha)
        self.check_alpha(alpha)
        try:
            df_alpha = self.alpha_dfs[alpha]
            dphi_ = df_alpha.dot(self.p)
        except KeyError:
            self.log(f"\tEvaluating gradient for alpha={alpha:.6f}")
            df_alpha = self.df(self.x0 + alpha * self.p)
            self.df_evals += 1
            self.alpha_dfs[alpha] = df_alpha
            dphi_ = df_alpha.dot(self.p)
            self.dphis[alpha] = dphi_
        return dphi_

    def got_alpha_phi_dphi(self, alpha):
        return (alpha in self.alpha_fs) and (alpha in self.alpha_dfs)

    def get_phi_dphi(self, what, alpha, check=True):
        """Wrapper that handles function/gradient evaluations."""
        alpha = float(alpha)
        whats = "f g fg".split()
        assert what in whats
        calc_funcs = {
            "f": self._phi,
            "g": self._dphi,
        }
        result = [calc_funcs[w](alpha) for w in what]
        # Check if we got both phi and dphi for alpha now. If so we
        # can check if the chosen condition (Wolfe/approx. Wolfe) is
        # satisfied.
        if (
            check
            and (alpha > 0.0)
            and self.can_eval_cond_func(alpha)
            and self.cond_func(alpha)
        ):
            self.log(f"Line search condition is satisfied for α={alpha:.6f}.")
            raise LineSearchConverged(alpha)
        # Dont return a list if only f or g was requested.
        if len(what) == 1:
            result = result[0]
        return result

    def get_fg(self, what, alpha):
        """Lookup raw function/gradient values for a given alpha."""
        whats = "f g fg".split()
        assert what in whats
        lookups = {
            "f": self.alpha_fs,
            "g": self.alpha_dfs,
        }
        result = [lookups[w][alpha] for w in what]
        if len(what) == 1:
            result = result[0]
        return result

    def sufficiently_decreased(self, alpha):
        """Sufficient decrease/Armijo condition."""
        return self._phi(alpha) <= (self.phi0 + self.c1 * alpha * self.dphi0)

    def curvature_condition(self, alpha):
        return self._dphi(alpha) >= self.c2 * self.dphi0

    def strong_curvature_condition(self, alpha):
        return abs(self._dphi(alpha)) <= -self.c2 * self.dphi0

    def wolfe_condition(self, alpha):
        """Normal, not strong, Wolfe condition."""
        return self.sufficiently_decreased(alpha) and self.curvature_condition(alpha)

    def strong_wolfe_condition(self, alpha):
        """Strong wolfe condition."""
        return self.sufficiently_decreased(alpha) and self.strong_curvature_condition(
            alpha
        )

    @abc.abstractmethod
    def run_line_search(self):
        raise NotImplementedError

    def run(self):
        self.prepare_line_search()

        # Normal termination
        try:
            self.log(f"Starting {self.__class__.__name__} line search.")
            alpha = self.run_line_search()
        # Termination in get_phi_dphi
        except LineSearchConverged as lsc:
            alpha = lsc.alpha
        # Failed LineSearch
        except LineSearchNotConverged:
            alpha = None

        result = LineSearchResult(
            converged=bool(alpha),
            alpha=alpha,
            # The gradient at the endpoint of the line search may not
            # have been evaluted, but the function value was always
            # evaluated, except when the line search did not converge.
            f_new=self.alpha_fs.get(alpha, None),
            g_new=self.alpha_dfs.get(alpha, None),
            f_evals=self.f_evals,
            df_evals=self.df_evals,
            dphi0=self.dphi0,
        )
        return result
