#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import fit_rigid, procrustes
from pysisyphus.optimizers.Optimizer import Optimizer

# [1] Nocedal, Wright - Numerical Optimization, 2006


class LBFGS(Optimizer):
    def __init__(self, geometry, alpha=1.0, keep_last=7, beta=1, **kwargs):
        self.alpha = alpha
        self.beta = beta
        assert isinstance(keep_last, int) and keep_last > 0
        self.keep_last = keep_last

        super().__init__(geometry, **kwargs)

        self.steps_ = list()
        self.grad_diffs = list()

    def restrict_step(self, steps):
        too_big = np.abs(steps) > self.max_step
        self.log(f"Found {np.sum(too_big)} big step components.")
        signs = np.sign(steps[too_big])
        # import pdb; pdb.set_trace()
        steps[too_big] = signs * self.max_step
        return steps

    def bfgs_multiply(self, s_list, y_list, force, beta=1):
        """Get a L-BFGS step.
        
        Algorithm 7.4 Nocedal, Num. Opt., p. 178."""
        q = -force
        cycles = len(s_list)
        alphas = list()
        rhos = list()
        # Store rho and alphas as they are also needed in the second loop
        for i in reversed(range(cycles)):
            s = s_list[i]
            y = y_list[i]
            rho = 1/y.dot(s)
            rhos.append(rho)
            alpha = rho * s.dot(q)
            alphas.append(alpha)
            q = q - alpha*y
        # Restore original order, so that rho[i] = 1/s_list[i].dot(y_list[i]) etc.
        alphas = alphas[::-1]
        rhos = rhos[::-1]

        if cycles > 0:
            s = s_list[-1]
            y = y_list[-1]
            gamma = s.dot(y) / y.dot(y)
            r = gamma * q
        else:
            r = beta * q

        for i in range(cycles):
            s = s_list[i]
            y = y_list[i]
            beta = rhos[i] * y.dot(r)
            r = r + s*(alphas[i] - beta)

        return r

    def optimize(self):
        if self.is_cos and self.align:
            rot_vecs, rot_vec_lists, _ = fit_rigid(
                self.geometry,
                vector_lists=(self.steps, self.forces, self.steps_, self.grad_diffs)
            )
            rot_steps, rot_forces, rot_steps_, rot_grad_diffs = rot_vec_lists
            self.steps = rot_steps
            self.forces = rot_forces
            self.steps_ = rot_steps_
            self.grad_diffs = rot_grad_diffs

        forces = self.geometry.forces
        self.forces.append(forces)
        energy = self.geometry.energy
        self.energies.append(energy)

        if self.cur_cycle > 0:
            prev_forces = self.forces[-2]
            grad_diff = prev_forces - forces
            self.grad_diffs.append(grad_diff)

        step = -self.bfgs_multiply(self.steps_, self.grad_diffs, forces, beta=self.beta)
        # step = self.scale_by_max_step(step)
        # norm = np.linalg.norm(step)
        # if norm > self.max_step:
            # step = step / norm * self.max_step
        step = self.restrict_step(step)
        # Only keep 'keep_last' cycles
        self.steps_ = self.steps.copy()[-self.keep_last:]
        self.grad_diffs = self.grad_diffs[-self.keep_last:]
        return step
