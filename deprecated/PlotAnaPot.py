import matplotlib.pyplot as plt
import numpy as np

class PlotAnaPot:

    def __init__(self, geometry, xlim, ylim, levels):
        self.geometry = geometry
        self.xlim = xlim
        self.ylim = ylim
        self.levels = levels

        self.fig, self.ax = plt.subplots(figsize=(8,8))

        x = np.linspace(*self.xlim, 100)
        y = np.linspace(*self.ylim, 100)
        X, Y = np.meshgrid(x, y)
        atoms = self.geometry.atoms
        pot_coords = np.stack((X, Y))
        pot = self.geometry.calculator.get_energy(atoms, pot_coords)["energy"]
        levels = np.linspace(*self.levels)
        contours = self.ax.contour(X, Y, pot, levels)

    def plot(self, coords):
        self.ax.plot(*zip(*coords), "ro", ls="-")

    def plot_gs(self, coords, pivot_coords, micro_coords):
        # Pivot points
        self.ax.plot(*zip(*pivot_coords), "bo", ls="-", label="pivot")
        # Constrained optmizations
        for mc in micro_coords:
            self.ax.plot(*zip(*mc), "yo", ls="-")
            for i, m in enumerate(mc):
                self.ax.text(*m, str(i))
        self.plot(coords)
        plt.legend()

    def show(self):
        plt.show()
