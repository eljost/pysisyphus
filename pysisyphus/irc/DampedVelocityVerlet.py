from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter

import matplotlib.pyplot as plt
import numpy as np

# [1] Following Reaction Pathways Using a Damped Classical Trajectory Algorithm
# 10.1021/jp012125b


class DampedVelocityVerlet(IRC):

    def __init__(self, geometry, v0=0.06, dt0=0.1, error_tol=0.003,
                 energy_lowering=1e-4, **kwargs):
        super(DampedVelocityVerlet, self).__init__(geometry, **kwargs)

        self.v0 = v0
        self.error_tol = error_tol
        self.dt0 = dt0

        self.mm_inv = self.geometry.mm_inv

        step_header = "damping dt dt_new error".split()
        step_fmts = [".2f", ".4f", ".4f", ".3E"]
        self.step_formatter = TableFormatter(step_header, step_fmts, 10)

    def prepare(self, direction):
        # In contrast to the paper we don't start from the TS geometry but
        # instead do an initial displacement along the imaginary mode.
        # We still consider the TS geometry and the transition vector.
        super(DampedVelocityVerlet, self).prepare(direction)

        #self.irc_coords = [self.ts_coords]
        init_factor = 1 if (direction == "forward") else -1
        initial_velocity, _ = self.damp_velocity(
                                        init_factor * self.transition_vector
        )
        self.velocities = [initial_velocity]
        acceleration = self.mm_inv.dot(self.geometry.forces)
        self.accelerations = [acceleration]
        self.time_steps = [self.dt0]

    def damp_velocity(self, velocity):
        # Eq. (4) in [1]
        damping_factor = self.v0 / np.linalg.norm(velocity)
        damped_velocity = velocity * damping_factor
        return damped_velocity, damping_factor

    def estimate_error(self):
        # See Fig. 1 and Eq. (5) in [1]
        time_step_sum = self.time_steps[-2] + self.time_steps[-1]
        # x'
        ref_coords = (
            self.irc_coords[-2]
            + self.velocities[-2]*time_step_sum
            + 0.5*self.accelerations[-2]*time_step_sum**2
        )
        coords_diff = np.abs(self.geometry.coords - ref_coords)
        largest_component = np.max(np.abs(coords_diff))
        coords_diff_norm = np.linalg.norm(coords_diff)
        estimated_error = max((largest_component, coords_diff_norm))
        return estimated_error

    def step(self):
        last_acceleration = self.accelerations[-1]
        last_velocity = self.velocities[-1]
        time_step = self.time_steps[-1]

        # Get new acceleration
        acceleration = self.mm_inv.dot(self.geometry.forces)
        self.accelerations.append(acceleration)
        self.irc_coords.append(self.geometry.coords)
        self.irc_energies.append(self.geometry.energy)
        # Calculate new coords and velocity
        # Eq. (2) in [1]
        coords = (self.geometry.coords
                  + last_velocity*time_step
                  + 0.5*last_acceleration*time_step**2
        )
        self.geometry.coords = coords
        self.irc_coords.append(coords)
        # Update velocity
        velocity = (last_velocity
                    + 0.5*(last_acceleration+acceleration)*time_step
        )
        # Damp velocity
        damped_velocity, damping_factor = self.damp_velocity(velocity)
        self.velocities.append(damped_velocity)

        if self.cur_step is 0:
            # No error estimated after the first step
            estimated_error = self.error_tol
        else:
            estimated_error = self.estimate_error()

        # Get next time step from error estimation
        new_time_step = time_step * np.abs(self.error_tol
                                           / estimated_error)**(1/3)
        # Constrain time step between 0.0025 fs <= time_step <= 3.0 fs
        new_time_step = min(new_time_step, 3)
        new_time_step = max(new_time_step, 0.025)
        self.time_steps.append(new_time_step)

        print(self.step_formatter.header)
        print(self.step_formatter.line(damping_factor, time_step,
                                       new_time_step, estimated_error)
        )
