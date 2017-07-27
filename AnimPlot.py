#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class AnimPlot:

    def __init__(self, calculator, optimizer,
                 xlim=(-1, 1), ylim=(-1, 1), num=100,
                 figsize=(8, 8), levels=(-150, 5, 30)):

        self.calculator = calculator
        self.optimizer = optimizer

        coords = self.optimizer.coords
        coords = [c.reshape(-1, 3) for c in coords]
        self.coords = coords
        forces = optimizer.forces
        forces = [f.reshape((-1, 3)) for f in forces]
        self.forces = forces

        self.fig, self.ax = plt.subplots(figsize=figsize)

        # Calculate the potential
        x = np.linspace(*xlim, 100)
        y = np.linspace(*ylim, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y, Z))
        pot = self.calculator.get_energy(fake_atoms, pot_coords)["energy"]

        # Draw the potentials contourlines
        levels = np.linspace(*levels)
        contours = self.ax.contour(X, Y, pot, levels)
        self.ax.clabel(contours, inline=1, fontsize=5)

        images_x = self.coords[0][:,0]
        images_y = self.coords[0][:,1]
        forces_x = self.forces[0][:,0]
        forces_y = self.forces[0][:,1]
        # Create artists, so we can update their data later
        self.images, = self.ax.plot(images_x, images_y, "ro", ls="-")
        self.quiv = self.ax.quiver(images_x, images_y, forces_x, forces_y)

    def func(self, frame):
        self.fig.suptitle("Cycle {}".format(frame))

        images_x = self.coords[frame][:,0]
        images_y = self.coords[frame][:,1]
        self.images.set_xdata(images_x)
        self.images.set_ydata(images_y)

        forces_x = self.forces[frame][:,0]
        forces_y = self.forces[frame][:,1]
        offsets = np.stack((images_x, images_y), axis=-1).flatten()
        # https://stackoverflow.com/questions/19329039
        # https://stackoverflow.com/questions/17758942
        self.quiv.set_offsets(offsets)
        self.quiv.set_UVC(forces_x, forces_y)

        return self.images, self.quiv

    def animate(self):
        cycles = range(self.optimizer.cur_cycle)
        an = animation.FuncAnimation(self.fig,
                                     self.func,
                                     frames=cycles,
                                     interval=500)
        plt.show()
