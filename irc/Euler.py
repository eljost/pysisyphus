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
        self.irc_energies.append(energy)
        energy_diff = self.irc_energies[0] - energy

        print(self.step_formatter.header)
        print(self.step_formatter.line(energy, energy_diff, gradient_norm))

        self.coords -= self.step_length*gradient/gradient_norm

    def show2d(self):
        fig, ax = plt.subplots(figsize=(8, 8))

        xlim = (-1.75, 1.25)
        ylim = (-0.5, 2.25)
        levels = (-150, -15, 40)
        x = np.linspace(*xlim, 100)
        y = np.linspace(*ylim, 100)
        X, Y = np.meshgrid(x, y)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y))
        pot = self.geometry.calculator.get_energy(fake_atoms,
                                                  pot_coords)["energy"]
        levels = np.linspace(*levels)
        contours = ax.contour(X, Y, pot, levels)

        ax.plot(*zip(*self.all_coords), "ro", ls="-")
        plt.show()
