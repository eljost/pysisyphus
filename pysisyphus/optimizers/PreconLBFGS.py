#!/usr/bin/env python3

from collections import deque

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.line_searches import backtracking
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.precon import get_precon


class PreconLBFGS(Optimizer):

    def __init__(self, geometry, history=7, **kwargs):
        assert geometry.coord_type == "cart"
        super().__init__(geometry, **kwargs)

        self.history = history

        self.grad_diffs = deque(maxlen=self.history)
        self.steps_ = deque(maxlen=self.history)

        int_ = RedundantCoords(self.geometry.atoms, self.geometry.cart_coords)
        self.bonds = int_.bond_indices
        # self.bends = int_.bending_indices
        self.bends = list()

    def get_precon(self):
        return get_precon(self.geometry.atoms, self.geometry.coords,
                          self.bonds, self.bends)

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        # import pdb; pdb.set_trace()
        if self.cur_cycle > 0:
            self.grad_diffs.append(-forces - -self.forces[-2])
            self.steps_.append(self.steps[-1])

            P = self.get_precon()
            # import pdb; pdb.set_trace()
            step = -bfgs_multiply(self.steps_, self.grad_diffs, forces, P=P)
        else:
            step = forces

        step_dir = step / np.linalg.norm(step)

        f = lambda coords: self.geometry.get_energy_at(coords)
        df = lambda coords: self.geometry.get_energy_and_forces_at(coords)["forces"]

        kwargs = {
            "f": f,
            "df": df,
            "x0": self.geometry.coords,
            "p": step_dir,
            "f0": energy,
            "g0": -forces,
            "alpha_init": 1.,
        }
        alpha, f_new, g_new = backtracking(**kwargs)
        step = alpha * step_dir
        return step
