from collections import deque
from functools import partial
from typing import Literal, Optional

import numpy as np
from scipy.sparse.linalg import spsolve

from pysisyphus.calculators import Dimer
from pysisyphus.cos.GrowingNT import GrowingNT
from pysisyphus.Geometry import Geometry
from pysisyphus.line_searches import *
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.precon import precon_getter, PreconKind


LineSearch = Literal["armijo", "armijo_fg", "strong_wolfe", "hz", None, False]


class PreconLBFGS(Optimizer):
    def __init__(
        self,
        geometry: Geometry,
        alpha_init: float = 1.0,
        history: int = 7,
        precon: bool = True,
        precon_update: int = 1,
        precon_getter_update: Optional[int] = None,
        precon_kind: PreconKind = "full",
        max_step_element: Optional[float] = None,
        line_search: LineSearch = "armijo",
        c_stab: Optional[float] = None,
        **kwargs,
    ) -> None:
        """Preconditioned limited-memory BFGS optimizer.

        See pysisyphus.optimizers.precon for related references.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        alpha_init
            Initial scaling factor for the first trial step in the excplicit line search.
        history
            History size. Keep last 'history' steps and gradient differences.
        precon
            Wheter to use preconditioning or not.
        precon_update
            Recalculate preconditioner P in every n-th cycle with the same topology.
        precon_getter_update
            Recalculate topology for preconditioner P in every n-th cycle. It is usually
            sufficient to only determine the topology once at the beginning.
        precon_kind
            What types of primitive internal coordinates to consider in the preconditioner.
        max_step_element
            Maximum component of the absolute step vector when no line search is carried out.
        line_search
            Whether to use explicit line searches and if so, which kind of line search.
        c_stab
            Regularization constant c in (H + cI)⁻¹ in atomic units.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Optimizer baseclass.
        """
        if precon:
            assert geometry.coord_type in (
                "cart",
                "cartesian",
            ), "Preconditioning makes only sense with 'coord_type: cart|cartesian'"
        # Set before calling the superclass constructor so we can use the value
        # in _set_opt_restart_info.
        self.history = history

        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init
        self.precon = precon
        # Disable preconditioning for 1 atom species, e.g., analytical potentials
        if self.precon and (self.geometry.cart_coords.size == 3):
            self.precon = None
        self.precon_update = precon_update
        self.precon_getter_update = precon_getter_update
        self.precon_kind = precon_kind

        self.P = None
        try:
            is_dimer = isinstance(self.geometry.calculator, Dimer)
        # COS objectes may not have a calculator
        except AttributeError:
            is_dimer = False
        go_uphill = is_dimer or isinstance(self.geometry, GrowingNT)

        if c_stab is None:
            self.log("No c_stab specified.")
            if go_uphill:
                c_stab = 0.103  # 1 eV/Å²
                self.log(
                    "Found climbing calculation. Using higher regularization "
                    f"c_stab={c_stab:.4f} au/bohr**2."
                )
            else:
                c_stab = 0.0103  # 0.1 eV/Å²

        self.c_stab = float(c_stab)
        self.log(f"Using c_stab={self.c_stab:.6f}")

        if max_step_element is None and go_uphill:
            max_step_element = 0.25
            self.log(
                f"Found climbing calculation. Using max_step_element={max_step_element:.2f}"
            )
        elif max_step_element is None and line_search in (None, False):
            max_step_element = 0.2
        self.max_step_element = max_step_element

        # Disable linesearch if max_step_element is set
        if self.max_step_element is not None:
            self.log(
                f"max_step_element={max_step_element:.6f} given. Setting "
                "line_search to 'None'."
            )
            line_search = None

        self.line_search = line_search
        ls_cls = {
            "armijo": Backtracking,
            "armijo_fg": partial(Backtracking, use_grad=True),
            "strong_wolfe": StrongWolfe,
            "hz": HagerZhang,
            None: None,
            False: None,
        }
        self.line_search_cls = ls_cls[self.line_search]

        if not self.restarted:
            self.grad_diffs = deque(maxlen=self.history)
            self.steps_ = deque(maxlen=self.history)
            self.f_evals = 0
            self.df_evals = 0

    def get_precon_getter(self):
        return precon_getter(
            self.geometry, c_stab=self.c_stab, kind=self.precon_kind, logger=self.logger
        )

    def prepare_opt(self):
        if self.precon:
            self.precon_getter = self.get_precon_getter()

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "c_stab": self.c_stab,
            "grad_diffs": np.array(self.grad_diffs).tolist(),
            "steps_": np.array(self.steps_).tolist(),
            "precon_kind": self.precon_kind,
            "f_evals": self.f_evals,
            "df_evals": self.df_evals,
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.grad_diffs = deque(
            [np.array(gd) for gd in opt_restart_info["grad_diffs"]], maxlen=self.history
        )
        self.steps_ = deque(
            [np.array(cd) for cd in opt_restart_info["steps_"]], maxlen=self.history
        )

        self.c_stab = opt_restart_info["c_stab"]
        self.precon_kind = opt_restart_info["precon_kind"]
        self.precon_getter = self.get_precon_getter()
        self.f_evals = opt_restart_info["f_evals"]
        self.df_evals = opt_restart_info["df_evals"]

    def scale_max_element(self, step, max_step_element):
        max_element = np.abs(step).max()
        if max_element > max_step_element:
            step *= max_step_element / max_element
        return step

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        norm = np.linalg.norm(forces)
        if not self.is_cos:
            fmt = " >12.6f"
            self.log(f" Current energy={energy:{fmt}} au")
            if self.cur_cycle > 0:
                prev_energy = self.energies[-2]
                self.log(f"Previous energy={prev_energy:{fmt}} au")
                self.log(f"             ΔE={energy-prev_energy:{fmt}} au")
        self.log(f"norm(forces)={norm:.6f}")

        # Check for preconditioner getter update
        if (
            self.precon_getter_update
            and self.cur_cycle > 0
            and self.cur_cycle % self.precon_getter_update == 0
        ):
            self.precon_getter = self.get_precon_getter()
            self.log("Calculated new preconditioner getter")

        # If requested, construct preconditioner
        if self.precon and (self.cur_cycle % self.precon_update == 0):
            self.P = self.precon_getter(self.geometry.coords)
            self.log("Calculated new preconditioner P")

        # Determine step direction

        # Preconditioned LBFGS in later cycles
        if self.cur_cycle > 0:
            self.grad_diffs.append(-forces - -self.forces[-2])
            self.steps_.append(self.steps[-1])
            step = bfgs_multiply(self.steps_, self.grad_diffs, forces, P=self.P)
        # Preconditioned steepest descent in cycle 0
        elif self.precon:
            step = spsolve(self.P, forces)
        # Steepest descent w/o preconditioner in cycle 0
        else:
            step = forces

        step_norm = np.linalg.norm(step)
        self.log(f"norm(unscaled step)={step_norm:.6f} au")

        # Determine step length
        if self.line_search_cls is not None:
            kwargs = {
                "geometry": self.geometry,
                "p": step,
                "f0": energy,
                "g0": -forces,
                "alpha_init": self.alpha_init,
            }
            line_search = self.line_search_cls(**kwargs)
            line_search_result = line_search.run()
            try:
                alpha = line_search_result.alpha
                step *= alpha
                f_evals = line_search_result.f_evals
                df_evals = line_search_result.df_evals
                self.log(
                    f"Evaluated {f_evals} energies and {df_evals} gradients in line search."
                )
                self.f_evals += f_evals
                self.df_evals += df_evals
            except TypeError:
                self.log("Line search did not converge!")
                step = self.scale_max_element(step, self.max_step_element)
        else:
            step = self.scale_max_element(step, self.max_step_element)
        step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={step_norm:.6f} au")
        return step
