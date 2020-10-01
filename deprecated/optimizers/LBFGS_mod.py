#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import fit_rigid, procrustes
from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer
from pysisyphus.optimizers.closures import bfgs_multiply

# [1] Nocedal, Wright - Numerical Optimization, 2006

import matplotlib.pyplot as plt

def plot_steps_hist(steps, cycle, title):
    fig, ax = plt.subplots()
    ax.hist(steps)
    fig.suptitle(f"Cycle {cycle:03d}, {title}")
    fn = f"cycle_{cycle:03d}_{title}.png"
    fig.savefig(fn)
    plt.close()


def check_step(geom, steps, cycle):
    imgs = geom.images
    nimgs = len(imgs)
    shape = (nimgs, -1, 3)
    coords = geom.coords.reshape(*shape)
    steps_ = steps.reshape(*shape)
    ind = np.abs(steps_).argmax()
    ind = np.unravel_index(ind, steps_.shape)
    max_img = imgs[ind[0]]
    fn = f"max_img_cycle_{cycle:03d}.xyz"
    with open(fn, "w") as handle:
        handle.write(max_img.as_xyz(ind))


class LBFGS(BacktrackingOptimizer):
    def __init__(self, geometry, alpha=1.0, keep_last=15, bt_force=20, **kwargs):
        self.keep_last = keep_last
        super().__init__(geometry, alpha=alpha, bt_force=bt_force, **kwargs)

        self.sigmas = list()
        self.grad_diffs = list()

    def prepare_opt(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)
        # Calculate initial forces before the first iteration
        self.coords.append(self.geometry.coords)
        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)

    def scale_by_max_step(self, steps):
        steps_max = np.abs(steps).max()
        step_norm = np.linalg.norm(steps)
        self.log(f"Unscaled norm(step)={step_norm:.4f}")
        # check_step(self.geometry, steps, self.cur_cycle)
        # plot_steps_hist(steps, self.cur_cycle, "0_Unscaled")
        if steps_max > self.max_step:
            fact = self.max_step / steps_max
            self.log(f"Scaling step with factor={fact:.4f}")
            steps *= self.max_step / steps_max
            step_norm = np.linalg.norm(steps)
            self.log(f"Scaled norm(step)={step_norm:.4f}")
            plot_steps_hist(steps, self.cur_cycle, "1_Scaled")
        return steps

    def optimize(self):
        prev_coords = self.coords[-1]
        prev_forces = self.forces[-1]

        step = bfgs_multiply(self.sigmas, self.grad_diffs, prev_forces)
        step = self.scale_by_max_step(step)

        new_coords = prev_coords + self.alpha*step

        coords_tmp = prev_coords.copy()
        forces_tmp = prev_forces.copy()

        self.geometry.coords = new_coords
        if self.is_cos and self.align:
            rot_vecs, rot_vec_lists, _ = fit_rigid(
                self.geometry,
                (prev_coords, prev_forces),
                vector_lists=(self.sigmas, self.grad_diffs)
            )
            prev_coords, prev_forces = rot_vecs
            # self.sigmas, self.grad_diffs = rot_vec_lists
            rot_sigmas, rot_grad_diffs = rot_vec_lists
            # if sigs:
                # import pdb; pdb.set_trace()
            np.testing.assert_allclose(np.linalg.norm(rot_sigmas),
                                       np.linalg.norm(self.sigmas)
            )
            np.testing.assert_allclose(np.linalg.norm(rot_grad_diffs),
                                       np.linalg.norm(self.grad_diffs)
            )
            self.sigmas = rot_sigmas
            self.grad_diffs = rot_grad_diffs

        new_forces = self.geometry.forces
        new_energy = self.geometry.energy

        skip = self.backtrack(new_forces, prev_forces)
        print("alpha", self.alpha)
        if skip:
            self.geometry.coords = coords_tmp
            return None

        sigma = new_coords - prev_coords
        self.sigmas.append(sigma)
        grad_diff = prev_forces - new_forces
        self.grad_diffs.append(grad_diff)

        # if len(self.sigmas) == self.keep_last:
             # import pdb; pdb.set_trace()
        self.sigmas = self.sigmas[-self.keep_last:]
        self.grad_diffs = self.grad_diffs[-self.keep_last:]

        # Because we add the step later on we restore the original
        # coordinates and set the appropriate energies and forces.
        self.geometry.coords = prev_coords
        self.geometry.forces = new_forces
        self.geometry.energy = new_energy

        self.forces.append(new_forces)
        self.energies.append(new_energy)

        return step

    # def save_also(self):
        # return {
            # "alpha": self.alpha,
        # }
