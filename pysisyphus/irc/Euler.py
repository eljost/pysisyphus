import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter


class Euler(IRC):

    def __init__(self, geometry, step_length=0.01, **kwargs):
        super(Euler, self).__init__(geometry, step_length, **kwargs)

        step_header = "E/au ΔE(TS)/au |gradient|".split()
        step_fmts = [".4f", ".4f", ".4f"]
        self.step_formatter = TableFormatter(step_header, step_fmts, 10)

    def step(self):
        gradient = self.gradient
        gradient_norm = np.linalg.norm(gradient)
        energy = self.energy
        self.irc_coords.append(self.geometry.coords)
        self.irc_energies.append(energy)
        energy_diff = self.irc_energies[0] - energy

        print(self.step_formatter.header)
        print(self.step_formatter.line(energy, energy_diff, gradient_norm))

        self.coords -= self.step_length*gradient/gradient_norm
