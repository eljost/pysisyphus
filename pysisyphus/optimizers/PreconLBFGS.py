#!/usr/bin/env python3

from collections import deque

import numpy as np
from scipy.sparse.linalg import spsolve

from pysisyphus.calculators import Dimer
from pysisyphus.line_searches import *
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.precon import precon_getter


class PreconLBFGS(Optimizer):

    def __init__(self, geometry, alpha_init=1., history=7, precon=True,
                 max_step_element=None, line_search="armijo",
                 c_stab=None, **kwargs):
        assert geometry.coord_type == "cart", \
            "Preconditioning makes only sense with 'coord_type: cart'"
        super().__init__(geometry, **kwargs)

        self.alpha_init = alpha_init
        self.history = history
        self.precon = precon

        is_dimer = isinstance(self.geometry.calculator, Dimer)
        if c_stab is None:
            self.log("No c_stab specified.")
            if is_dimer:
                self.log( "Found Dimer calculator, using bigger stabilization.")
                c_stab = 0.103  # 1 eV/Å²
            else:
                c_stab = 0.0103 # 0.1 eV/Å²

        self.c_stab = float(c_stab)
        self.log("Using c_stab={self.c_stab:.6f}")

        if max_step_element is None and is_dimer:
            max_step_element = 0.25
            self.log("Found Dimer calculator. Using "
                     "max_step_element={max_step_element:.2f}")
        self.max_step_element = max_step_element
        if self.max_step_element is not None:
            self.log(f"max_step_element={max_step_element:.6f} given. Setting "
                      "line_search to 'None'.")
            line_search = None

        self.line_search = line_search
        ls_cls = {
            "armijo": Backtracking,
            "strong_wolfe": StrongWolfe,
            "hz": HagerZhang,
            None: None,
            False: None,
        }
        self.line_search_cls = ls_cls[self.line_search]

        self.grad_diffs = deque(maxlen=self.history)
        self.steps_ = deque(maxlen=self.history)

    def prepare_opt(self):
        if self.precon:
            self.precon_getter = precon_getter(self.geometry, c_stab=self.c_stab)

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        # Steepest descent fallback
        step = forces

        # Construct preconditoner if requested
        P = None
        if self.precon:
            P = self.precon_getter(self.geometry.coords)
            step = spsolve(P, forces)

        if self.cur_cycle > 0:
            self.grad_diffs.append(-forces - -self.forces[-2])
            self.steps_.append(self.steps[-1])
            step = -bfgs_multiply(self.steps_, self.grad_diffs, forces, P=P)

        step_dir = step / np.linalg.norm(step)

        if self.line_search_cls is not None:
            kwargs = {
                "geometry": self.geometry,
                "p": step_dir,
                "f0": energy,
                "g0": -forces,
                "alpha_init": self.alpha_init,
            }
            line_search = self.line_search_cls(**kwargs)
            line_search_result = line_search.run()
            alpha = line_search_result.alpha
            step = alpha * step_dir
        else:
            max_element = np.abs(step).max()
            if max_element > self.max_step_element:
                step *= self.max_step_element / max_element
        return step
