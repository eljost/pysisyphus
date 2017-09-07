#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

class IRC:

    def __init__(self, geometry, max_step=0.1, max_cycles=50):
        self.geometry = geometry
        self.max_step = max_step
        self.max_cycles = max_cycles

        self.hessian = self.geometry.hessian

        assert(max_step > 0), "max_step has to be > 0"

        self.coords_list = list()

    # set hessian

    def run(self):
        last_energy = None
        i = 0
        self.energies = list()
        self.coords = [self.geometry.coords, ]
        while True:
            print(f"Step {i}")
            self.step()
            self.coords.append(self.geometry.coords)
            this_energy = self.geometry.energy
            self.energies.append(this_energy)
            if i == self.max_cycles:
                break
            elif last_energy and (this_energy > last_energy):
                print("Energy increased!")
                break
            elif last_energy and abs(last_energy - this_energy) <= 1e-4:
                print("Energy converged!")
                break
            last_energy = this_energy
            i += 1

        self.coords = np.array(self.coords)
        self.postprocess()

    def postprocess(self):
        pass
