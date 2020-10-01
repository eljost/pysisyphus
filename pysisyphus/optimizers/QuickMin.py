import numpy as np

from pysisyphus.helpers import fit_rigid
from pysisyphus.optimizers.Optimizer import Optimizer


class QuickMin(Optimizer):

    def __init__(self, geometry, dt=0.35, **kwargs):
        super(QuickMin, self).__init__(geometry, **kwargs)

        self.dt = dt

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "velocity": self.velocities[-1].tolist(),
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        velocity = np.array(opt_restart_info["velocity"], dtype=float)
        self.velocities = [velocity, ]

    def prepare_opt(self):
        self.velocities = [np.zeros_like(self.geometry.coords), ]

    def reset(self):
        if self.coords[-1].size != self.coords[-2].size:
            self.log("Number of coordinates changed. Adapting velocity vector "
                     "to new coordinate number.")
            self.prepare_opt()
        self.log("Resetted optimizer")

    def optimize(self):
        if self.align and self.is_cos:
            (self.velocities[-1], ), _, _ = fit_rigid(self.geometry,
                                                      (self.velocities[-1], ))

        prev_velocities = self.velocities[-1]
        cur_forces = self.geometry.forces
        self.forces.append(cur_forces)
        self.energies.append(self.geometry.energy)

        norm = np.linalg.norm(cur_forces)
        if not self.is_cos:
            self.log(f"Current energy={self.energies[-1]:.6f}")
        self.log(f"norm(forces)={norm:.6f}")

        if self.cur_cycle == 0:
            tmp_velocities = np.zeros_like(prev_velocities)
        else:
            overlap = prev_velocities.dot(cur_forces)
            self.log(f"Overlap of previous and current forces: {overlap:.6f}")
            if overlap > 0:
                tmp_velocities = (overlap * cur_forces
                                  / cur_forces.dot(cur_forces))
            else:
                tmp_velocities = np.zeros_like(prev_velocities)
                self.log("Zeroed velocities")

        accelerations = cur_forces / self.geometry.masses_rep
        cur_velocities = tmp_velocities + self.dt*accelerations
        steps = cur_velocities*self.dt + 1/2*accelerations*self.dt**2
        steps = self.scale_by_max_step(steps)
        self.velocities.append(cur_velocities)
        velo_norm = np.linalg.norm(cur_velocities)
        acc_norm = np.linalg.norm(accelerations)
        self.log(f"norm(v) = {velo_norm:.4f}, norm(a) = {acc_norm:.4f}")

        return steps
