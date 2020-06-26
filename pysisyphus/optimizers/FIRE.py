import numpy as np

from pysisyphus.helpers import fit_rigid
from pysisyphus.optimizers.Optimizer import Optimizer

class FIRE(Optimizer):
    # https://doi.org/10.1103/PhysRevLett.97.170201

    def __init__(self, geometry, dt=0.1, dt_max=1, N_acc=2, f_inc=1.1,
                 f_acc=0.99, f_dec=0.5, n_reset=0, a_start=0.1, **kwargs):
        self.dt = dt
        self.dt_max = dt_max
        # Accelerate after N_acc cycles
        self.N_acc = N_acc
        self.f_inc = f_inc
        self.f_acc = f_acc
        self.f_dec = f_dec
        self.n_reset = n_reset
        self.a_start = a_start

        self.a = self.a_start
        # The current velocity
        self.v = np.zeros_like(geometry.coords)
        # Store the velocities for every step
        self.velocities = [self.v, ]
        self.time_deltas = [self.dt, ]

        super().__init__(geometry, **kwargs)

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "a": self.a,
            "dt": self.dt,
            "n_reset": self.n_reset,
            "time_delta": self.time_deltas[-1],
            "velocity": self.velocities[-1].tolist(),
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.a = opt_restart_info["a"]
        self.dt = opt_restart_info["dt"]
        self.n_reset = opt_restart_info["n_reset"]
        self.time_deltas = [opt_restart_info["time_delta"], ]
        velocity = np.array(opt_restart_info["velocity"], dtype=float)
        self.velocities = [velocity, ]

    def reset(self):
        pass

    def optimize(self):
        if self.is_cos and self.align:
            (self.v, ), _, _  = fit_rigid(self.geometry, (self.v, ))

        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)
        forces = self.forces[-1]
        mixed_v = (
            # As 'a' gets bigger we keep less old v.
            (1.0 - self.a) * self.v +
            # As 'a' gets bigger we use more new v
            # from the current forces.
             self.a * np.sqrt(
                        np.dot(self.v, self.v) / np.dot(forces, forces)
                      ) * forces
        )
        last_v = self.velocities[-1]
        # Check if forces are still aligned with the last velocity
        if (self.cur_cycle > 0) and (np.dot(last_v, forces) > 0):
            if self.n_reset > self.N_acc:
                self.dt = min(self.dt * self.f_inc, self.dt_max)
                self.a *= self.f_acc
            self.n_reset += 1
        else:
            # Reset everything when 'forces' and 'last_v' aren't
            # aligned anymore.
            mixed_v = np.zeros_like(forces)
            self.log("resetted velocities")
            self.a = self.a_start
            self.dt *= self.f_acc
            self.n_reset = 0

        v = mixed_v + self.dt * forces
        self.velocities.append(v)
        self.time_deltas.append(self.dt)
        steps = self.dt * v
        steps = self.scale_by_max_step(steps)

        velo_norm = np.linalg.norm(v)
        self.log(f"dt = {self.dt:.4f}, norm(v) {velo_norm:.4f}")

        return steps

    def __str__(self):
        return "FIRE optimizer"
