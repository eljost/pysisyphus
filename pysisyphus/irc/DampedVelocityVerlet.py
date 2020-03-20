import numpy as np

from pysisyphus.constants import BOHR2M, AU2J, AMU2KG
from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter


# [1] https://pubs.acs.org/doi/10.1021/jp012125b
#     Following Reaction Pathways Using a Damped Classical Trajectory Algorithm
#     Hratchian, Schlegel, 2002


class DampedVelocityVerlet(IRC):

    def __init__(self, geometry, v0=0.04, dt0=0.5, error_tol=0.003,
                 max_cycles=150, **kwargs):
        super().__init__(geometry, max_cycles=max_cycles, **kwargs)

        self.v0 = v0
        self.error_tol = error_tol
        self.dt0 = dt0

        step_header = "damping dt dt_new error".split()
        step_fmts = [".2f", ".4f", ".4f", ".3E"]
        self.step_formatter = TableFormatter(step_header, step_fmts, 10)

    def mw_grad_to_acc(self, mw_grad):
        """Takes care of the units for the mass-weighted gradient.
        Converts units of a mass-weighted gradient [Hartree/(Bohr*amu)]
        to units of acceleration [sqrt(amu)*Bohr/fs²]. The 1e30 comes
        from converting second² to femto second²."""
        return mw_grad * AU2J / AMU2KG / BOHR2M**2 / 1e30

    def prepare(self, direction):
        # In contrast to the paper we don't start from the TS geometry but
        # instead do an initial displacement along the imaginary mode as for
        # the other IRC classes. So there should always be a non-vanishing
        # gradient after the prepare call.
        super().prepare(direction)

        acceleration = -self.mw_grad_to_acc(self.mw_gradient)
        # init_velocity = acceleration
        init_velocity, _ = self.damp_velocity(acceleration)

        self.velocities = [init_velocity]
        self.accelerations = [acceleration]
        self.time_steps = [self.dt0]

    def damp_velocity(self, velocity):
        # Eq. (4) in [1]
        damping_factor = self.v0 / np.linalg.norm(velocity)
        self.log(f"Damping factor={damping_factor:.6f}")
        damped_velocity = velocity * damping_factor
        return damped_velocity, damping_factor

    def estimate_error(self, new_mw_coords):
        self.log("Error estimation")
        # See Fig. 1 and Eq. (5) in [1]
        cur_time_step = self.time_steps[-1]
        prev_time_step = self.time_steps[-2]
        time_step_sum = prev_time_step + cur_time_step
        self.log(f"\tSum of cur. and prev. timestep: {time_step_sum:.6f} fs")
        # Calculation of x' coordinates
        # irc_mw_coords
        #   -1: current coords  <- x_i-1 in the paper
        #   -2: previous coords <- x_i-2 in the paper
        # velocities/accelerations
        #   -1: new velo./acc. for next time step
        #   -2: velo./acc. that yielded the proposed coords x_i
        #   -3: velo./acc. that yielded previous structure <- v_i-2 and a_i-2 in the paper
        ref_coords = (
            self.irc_mw_coords[-2]
            # TODO: This should probably be -3 but it seems to give inferior results
            + self.velocities[-2]*time_step_sum
            + 0.5*self.accelerations[-2]*time_step_sum**2
        )
        diff = new_mw_coords - ref_coords
        diff /= np.sqrt(self.geometry.masses_rep)
        largest_component = np.max(np.abs(diff))
        norm = np.linalg.norm(diff)
        # Take either the largest component of the difference vector
        # or the norm of the difference vector.
        self.log("\tdiff=Proposed coords - Reference coords")
        self.log(f"\tmax(|diff|)={largest_component:.6f}")
        self.log(f"\ttnorm(diff)={norm:.6f}")
        # estimated_error = max(largest_component, norm)
        estimated_error = max(largest_component, norm)
        self.log(f"\testimated error={estimated_error:.6f}")
        return estimated_error

    def step(self):
        prev_acceleration = self.accelerations[-1]
        prev_velocity = self.velocities[-1]
        time_step = self.time_steps[-1]

        # Get new acceleration
        acceleration = -self.mw_grad_to_acc(self.mw_gradient)

        acc_normed = acceleration/np.linalg.norm(acceleration)
        prev_vel_normed = prev_velocity / np.linalg.norm(prev_velocity)
        ovlp = acc_normed @ prev_vel_normed
        self.log(f"a @ v_i-1={ovlp:.8f}")
        self.accelerations.append(acceleration)

        # Calculate new coords and velocity
        # Eq. (2) in [1]
        new_mw_coords = (self.mw_coords
                         + prev_velocity*time_step
                         + 0.5*prev_acceleration*time_step**2
        )

        # Update velocity
        velocity = prev_velocity + 0.5*(prev_acceleration+acceleration)*time_step
        # Damp velocity
        damped_velocity, damping_factor = self.damp_velocity(velocity)
        self.velocities.append(damped_velocity)

        if self.cur_cycle == 0:
            # No error estimation in the first step
            estimated_error = self.error_tol
        else:
            estimated_error = self.estimate_error(new_mw_coords)

        # Get next time step from error estimation
        new_time_step = time_step * (self.error_tol / estimated_error)**(1/3)
        # Constrain time step between 0.0025 fs and 3.0 fs
        new_time_step = min(new_time_step, 3.)
        new_time_step = max(new_time_step, 0.025)
        self.time_steps.append(new_time_step)
        self.log(f"\tCur. time step={time_step:.6f} fs")
        self.log(f"\tNext time step={new_time_step:.6f} fs")

        # Set new coordinates
        self.mw_coords = new_mw_coords
