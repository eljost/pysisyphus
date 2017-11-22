#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np


def get_coords_diffs(coords):
    cds = [0, ]
    for i in range(len(coords)-1):
        diff = np.linalg.norm(coords[i+1]-coords[i])
        cds.append(diff)
    cds = np.cumsum(cds)
    cds /= cds.max()
    return cds


class AnimPlot:

    def __init__(self, calculator, optimizer,
                 xlim=(-1, 1), ylim=(-1, 1), num=100,
                 figsize=(8, 8), levels=(-150, 5, 30),
                 interval=250):

        self.calculator = calculator
        self.optimizer = optimizer
        self.interval = interval

        self.coords = [c.reshape(-1, 3) for c in self.optimizer.coords]
        self.forces = [f.reshape((-1, 3)) for f in self.optimizer.forces]
        self.energies = self.optimizer.energies
        self.tangents = self.optimizer.tangents

        # ax: the contour plot
        # ax1: energy along the path
        self.fig, (self.ax, self.ax1) = plt.subplots(2, figsize=figsize,
                                        gridspec_kw = {'height_ratios':[3, 1]})

        self.pause = True
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

        # Calculate the potential
        x = np.linspace(*xlim, 100)
        y = np.linspace(*ylim, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y, Z))
        pot = self.calculator.get_energy(fake_atoms, pot_coords)["energy"]

        # Draw the contourlines of the potential
        levels = np.linspace(*levels)
        contours = self.ax.contour(X, Y, pot, levels)
        #self.ax.clabel(contours, inline=1, fontsize=5)
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")

        # Create a colorbar
        self.fig.subplots_adjust(right=0.8)
        cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
        self.fig.colorbar(contours, cax=cbar_ax)

        images_x = self.coords[0][:,0]
        images_y = self.coords[0][:,1]
        forces_x = self.forces[0][:,0]
        forces_y = self.forces[0][:,1]
        tangents_x = self.tangents[0][:,0]
        tangents_y = self.tangents[0][:,1]
        energies = self.energies[0]

        # Create artists, so we can update their data later
        # Image positions
        self.images, = self.ax.plot(images_x, images_y, "ro", ls="-")
        # Total forces
        self.total_forces_quiv = self.ax.quiver(images_x, images_y,
                                                forces_x, forces_y)
        # Tangents
        self.tangent_quiv = self.ax.quiver(images_x, images_y,
                                           tangents_x, tangents_y, color="b")

        # Energy along the path
        self.energies_plot, = self.ax1.plot(
            get_coords_diffs(self.coords[0]), energies, "ro", ls="-"
        )
        self.ax1.set_xlabel("q(x, y)")
        self.ax1.set_ylabel("f(x, y)")

    def func(self, frame):
        self.fig.suptitle("Cycle {}".format(frame))

        images_x = self.coords[frame][:,0]
        images_y = self.coords[frame][:,1]
        self.images.set_xdata(images_x)
        self.images.set_ydata(images_y)

        # Update total forces quiver
        forces_x = self.forces[frame][:,0]
        forces_y = self.forces[frame][:,1]
        offsets = np.stack((images_x, images_y), axis=-1).flatten()
        # https://stackoverflow.com/questions/19329039
        # https://stackoverflow.com/questions/17758942
        self.total_forces_quiv.set_offsets(offsets)
        self.total_forces_quiv.set_UVC(forces_x, forces_y)

        # Update tangent quiver
        tangents_x = self.tangents[frame][:,0]
        tangents_y = self.tangents[frame][:,1]
        self.tangent_quiv.set_offsets(offsets)
        self.tangent_quiv.set_UVC(tangents_x, tangents_y)

        coords_diffs = get_coords_diffs(self.coords[frame])
        energies = self.energies[frame]
        self.energies_plot.set_xdata(coords_diffs)
        self.energies_plot.set_ydata(energies)
        self.ax1.relim()
        self.ax1.autoscale_view()

    def animate(self):
        cycles = range(self.optimizer.cur_cycle)
        self.animation = animation.FuncAnimation(self.fig,
                                                 self.func,
                                                 frames=cycles,
                                                 interval=self.interval)
        plt.show()

    def on_keypress(self, event):
        """Pause on SPACE press."""
        #https://stackoverflow.com/questions/41557578
        if event.key == " ":
            if self.pause:
                self.animation.event_source.stop()
            else:
                self.animation.event_source.start()
            self.pause = not self.pause
