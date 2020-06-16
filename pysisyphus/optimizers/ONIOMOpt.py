import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class ONIOMOpt(Optimizer):

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        step = forces
        norm = np.linalg.norm(step)
        if norm > 0.5:
            direction = forces / np.linalg.forces(direction)
            step = 0.5 * direction
        return step
