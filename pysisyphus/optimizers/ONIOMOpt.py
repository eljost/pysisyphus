import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class ONIOMOpt(Optimizer):

    def __init__(self, geometry, *args, max_micro_cycles=None, **kwargs):
        super().__init__(geometry, *args, **kwargs)

        layers = self.geometry.layers
        print(f"found {len(layers)} layers: {layers}")

        if max_micro_cycles is None:
            max_micro_cycles = np.ones(len(layers), dtype=int)
            max_micro_cycles[0] = 5
        self.max_micro_cycles = max_micro_cycles
        print("max_micro_cycles", self.max_micro_cycles)

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        step = forces
        norm = np.linalg.norm(step)
        if norm > 0.5:
            direction = forces / np.linalg.norm(forces)
            step = 0.5 * direction
        return step
