#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.helpers import fit_rigid, procrustes
from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

# [1] Nocedal, Wright - Numerical Optimization, 2006

class BFGS(BacktrackingOptimizer):

    def __init__(self, geometry, alpha=1.0, bt_force=20, **kwargs):
        super(BFGS, self).__init__(geometry, alpha=alpha,
                                   bt_force=bt_force,
                                   **kwargs)

        self.eye = np.eye(len(self.geometry.coords))
        try:
            self.inv_hessian = self.geometry.get_initial_hessian()
        # ChainOfStates objects may not have get_initial_hessian
        except AttributeError:
            self.inv_hessian = self.eye.copy()
        if (hasattr(self.geometry, "internal")
            and (self.geometry.internal is not None)):
            raise Exception("Have to add hessian projections etc.")
        self.log("BFGS with align=True is somewhat broken right now, so "
                 "the images will be aligned only in the first iteration. "
        )

    def reset_hessian(self):
        self.inv_hessian = self.eye.copy()
        self.log("Resetted hessian")

    def prepare_opt(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)
        # Calculate initial forces before the first iteration
        self.coords.append(self.geometry.coords)
        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)

    def scale_by_max_step(self, steps):
        steps_max = np.abs(steps).max()
        if steps_max > self.max_step:
            fact = self.max_step / steps_max
            """
            fig, ax = plt.subplots()
            ax.hist(steps, bins=20)#"auto")
            title = f"max(steps)={steps_max:.04f}, fact={fact:.06f}"
            ax.set_title(title)
            l1 = ax.axvline(x=self.max_step, c="k")
            l2 = ax.axvline(x=-self.max_step, c="k")
            ax.add_artist(l1)
            ax.add_artist(l2)
            fig.savefig(f"cycle_{self.cur_cycle:02d}.png")
            plt.close(fig)
            """
            steps *= self.max_step / steps_max
        return steps

    def optimize(self):
        last_coords = self.coords[-1]
        last_forces = self.forces[-1]
        last_energy = self.energies[-1]

        unscaled_steps = self.inv_hessian.dot(last_forces)
        steps = self.scale_by_max_step(self.alpha*unscaled_steps)

        new_coords = last_coords + steps
        self.geometry.coords = new_coords

        # Hessian rotation seems faulty right now ...
        #if self.is_cos and self.align:
        #    (last_coords, last_forces, steps), _, self.inv_hessian = fit_rigid(
        #                                                    self.geometry,
        #                                                    (last_coords,
        #                                                     last_forces,
        #                                                     steps),
        #                                                    hessian=self.inv_hessian)

        new_forces = self.geometry.forces
        new_energy = self.geometry.energy
        skip = self.backtrack(new_forces, last_forces, reset_hessian=True)
        if skip:
            self.reset_hessian()
            self.geometry.coords = last_coords
            #self.scale_alpha(unscaled_steps, self.alpha)
            return None

        # Because we add the step later on we restore the original
        # coordinates and set the appropriate energies and forces.
        self.geometry.coords = last_coords
        self.geometry.forces = new_forces
        self.geometry.energy = new_energy

        self.forces.append(new_forces)
        self.energies.append(new_energy)
        # [1] Eq. 6.5, gradient difference, minus force difference
        y = -(new_forces - last_forces)
        sigma = new_coords - last_coords
        # [1] Eq. 6.7, curvature condition
        curv_cond = sigma.dot(y)
        if curv_cond < 0:
            self.log(f"curvature condition {curv_cond:.07} < 0!")
        rho = 1.0 / y.dot(sigma)
        if ((np.array_equal(self.inv_hessian, self.eye))
            # When align = True the above expression will evaluate to
            # False. So we also check if we are in the first iteration.
            or (self.cur_cycle == 0)):
            # [1] Eq. 6.20, p. 143
            beta = y.dot(sigma)/y.dot(y)
            self.inv_hessian = self.eye*beta
            self.log(f"Using initial guess for inverse hessian, beta={beta}")
        # Inverse hessian update
        A = self.eye - np.outer(sigma, y) * rho
        B = self.eye - np.outer(y, sigma) * rho
        self.inv_hessian = (A.dot(self.inv_hessian).dot(B)
                            + np.outer(sigma, sigma) * rho)

        return steps
