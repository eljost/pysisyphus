from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter

import matplotlib.pyplot as plt
import numpy as np

# [1] Following Reaction Pathways Using a Damped Classical Trajectory Algorithm
# 10.1021/jp012125b

class DampedVelocityVerlet(IRC):

    def __init__(self, geometry, v0=0.06, error_tol=0.003, **kwargs):
        super(DampedVelocityVerlet, self).__init__(geometry, **kwargs)

        self.v0 = v0
        self.error_tol = error_tol

        self.mass_mat_inv = np.linalg.inv(self.geometry.mass_mat)
        self.dvv_coords = [self.geometry.coords, self.geometry.coords]
        zeros = np.zeros_like(self.geometry.coords)
        import logging
        logging.warning("use transition vector for intial velocity")
        self.velocities = [zeros, zeros]
        acceleration = self.mass_mat_inv.dot(self.geometry.forces)
        #self.accelerations = [zeros, zeros]
        self.accelerations = [zeros, acceleration]
        self.time_steps = [0.1]

        step_header = "# damping dt dt_new error".split()
        step_fmts = ["d", ".2f", ".4f", ".4f", ".3E"]
        self.step_formatter = TableFormatter(step_header, step_fmts, 10)

    def step(self):
        last_acceleration = self.accelerations[-1]
        self.irc_energies.append(self.geometry.energy)
        last_velocity = self.velocities[-1]
        time_step = self.time_steps[-1]

        # Get new acceleration
        acceleration = self.mass_mat_inv.dot(self.geometry.forces)
        self.accelerations.append(acceleration)
        # Calculate new coords and velocity
        # Eq. (2) in [1]
        print("old coords", self.geometry.coords)
        new_coords = (self.geometry.coords
                      + last_velocity*time_step
                      + 0.5*last_acceleration*time_step**2
        )
        print("new coords", new_coords)
        # Update velocity
        new_velocity = (last_velocity
                        + 0.5*(last_acceleration+acceleration)*time_step
        )

        # Velocity damping
        # Eq. (4) in [1]
        damping_factor = self.v0 / np.linalg.norm(new_velocity)
        new_velocity *= damping_factor
        self.velocities.append(new_velocity)

        # Get next time step from error estimation
        # See Fig. 1 and Eq. (5) in [1]
        ref_coords = (
            self.dvv_coords[-2]
            + self.velocities[-2]*(self.time_steps[-1]+time_step)
            + 0.5*self.accelerations[-2]*(self.time_steps[-1]+time_step)**2
        )
        coords_diff = np.abs(new_coords - ref_coords)
        #print("coords_diff", coords_diff)
        largest_component = np.max(np.abs(coords_diff))
        #print("larg comp", largest_component)
        coords_diff_norm = np.linalg.norm(coords_diff)
        #print("coords diff norm", coords_diff_norm)
        estimated_error = max((largest_component, coords_diff_norm))
        new_time_step = time_step * np.linalg.norm(self.error_tol/estimated_error)**(1/3)
        # Constrain time step between 0.0025 fs <= time_step <= 3.0 fs
        new_time_step = min(new_time_step, 3)
        new_time_step = max(new_time_step, 0.025)
        self.time_steps.append(new_time_step)
        self.geometry.coords = new_coords

        print(self.step_formatter.header)
        print(self.step_formatter.line(self.cur_step, damping_factor,
                                    time_step, new_time_step, estimated_error))

    def show2d(self):
        fig, ax = plt.subplots(figsize=(8,8))

        xlim = (-1.75, 1.25)
        ylim = (-0.5, 2.25)
        levels=(-150, -15, 40)
        x = np.linspace(*xlim, 100)
        y = np.linspace(*ylim, 100)
        X, Y = np.meshgrid(x, y)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y))
        pot = self.geometry.calculator.get_energy(fake_atoms, pot_coords)["energy"]
        levels = np.linspace(*levels)
        contours = ax.contour(X, Y, pot, levels)

        ax.plot(*zip(*self.coords), "ro", ls="-")
        plt.legend()
        plt.show()
