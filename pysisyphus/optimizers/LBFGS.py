import numpy as np

from pysisyphus.helpers import fit_rigid
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.hessian_updates import double_damp
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.restrict_step import scale_by_max_step


# See pysisyphus.optimizers.hessian_updates for double damping ref.


class LBFGS(Optimizer):
    def __init__(self, geometry, keep_last=7, beta=1, max_step=0.2,
                 double_damp=True, gamma_mult=False, **kwargs):
        """[1] Nocedal, Wright - Numerical Optimization, 2006"""
        self.beta = beta
        self.keep_last = int(keep_last)
        self.double_damp = double_damp
        self.gamma_mult = gamma_mult

        self.coord_diffs = list()
        self.grad_diffs = list()

        super().__init__(geometry, max_step=max_step, **kwargs)

    def reset(self):
        self.coord_diffs = list()
        self.grad_diffs = list()

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "coord_diffs": np.array(self.coord_diffs).tolist(),
            "grad_diffs": np.array(self.grad_diffs).tolist(),
            "double_damp": self.double_damp,
            "gamma_mult": self.gamma_mult,
            "keep_last": self.keep_last,
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.coord_diffs = [np.array(cd) for cd in opt_restart_info["coord_diffs"]]
        self.grad_diffs = [np.array(gd) for gd in opt_restart_info["grad_diffs"]]
        for attr in ("double_damp", "gamma_mult", "keep_last"):
            setattr(self, attr, opt_restart_info[attr])

    def optimize(self):
        if self.is_cos and self.align:
            rot_vecs, rot_vec_lists, _ = fit_rigid(
                self.geometry,
                vector_lists=(self.steps, self.forces, self.coord_diffs, self.grad_diffs)
            )
            rot_steps, rot_forces, rot_coord_diffs, rot_grad_diffs = rot_vec_lists
            self.steps = rot_steps
            self.forces = rot_forces
            self.coord_diffs = rot_coord_diffs
            self.grad_diffs = rot_grad_diffs

        forces = self.geometry.forces
        self.forces.append(forces)
        energy = self.geometry.energy
        self.energies.append(energy)
        norm = np.linalg.norm(forces)
        self.log(f"norm(forces)={norm:.6f}")

        if self.cur_cycle > 0:
            y = self.forces[-2] - forces
            s = self.steps[-1]
            if self.double_damp:
                s, y = double_damp(s, y, s_list=self.coord_diffs, y_list=self.grad_diffs)
            self.grad_diffs.append(y)
            self.coord_diffs.append(s)

            # Drop superfluous oldest vectors
            self.coord_diffs = self.coord_diffs[-self.keep_last:]
            self.grad_diffs = self.grad_diffs[-self.keep_last:]

        step = bfgs_multiply(self.coord_diffs, self.grad_diffs, forces, beta=self.beta,
                             gamma_mult=self.gamma_mult, logger=self.logger)
        step = scale_by_max_step(step, self.max_step)

        return step
