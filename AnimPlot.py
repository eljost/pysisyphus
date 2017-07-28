#!/usr/bin/env python3


import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def get_coords_diffs(coords):
    cds = [np.linalg.norm(coords[0]-coords[i])
           for i in range(len(coords))
    ]
    cds = np.array(cds)
    cds /= cds.max()
    return cds

class AnimPlot:

    def __init__(self, calculator, optimizer,
                 xlim=(-1, 1), ylim=(-1, 1), num=100,
                 figsize=(8, 8), levels=(-150, 5, 30)):

        self.calculator = calculator
        self.optimizer = optimizer

        coords = self.optimizer.coords
        coords = [c.reshape(-1, 3) for c in coords]
        self.coords = coords
        self.energies = self.optimizer.energies
        forces = optimizer.forces
        forces = [f.reshape((-1, 3)) for f in forces]
        self.forces = forces

        self.fig, (self.ax, self.ax1) = plt.subplots(2, figsize=figsize,
                                        gridspec_kw = {'height_ratios':[3, 1]})
        self.pause = True
        self.fig.canvas.mpl_connect('key_press_event', self.on_click)

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
        #self.ax.clabel(contours, inline=1, fontsize=5)
        # How do you add a colorbar via the axis object?
        #plt.colorbar(contours)
        self.fig.subplots_adjust(right=0.8)
        cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
        self.fig.colorbar(contours, cax=cbar_ax)

        images_x = self.coords[0][:,0]
        images_y = self.coords[0][:,1]
        forces_x = self.forces[0][:,0]
        forces_y = self.forces[0][:,1]
        energies = self.energies[0]
        # Create artists, so we can update their data later
        self.images, = self.ax.plot(images_x, images_y, "ro", ls="-")
        self.quiv = self.ax.quiver(images_x, images_y, forces_x, forces_y)
        self.energies_plot, = self.ax1.plot(
            get_coords_diffs(self.coords[0]), energies, "ro", ls="-"
        )

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

        coords_diffs = get_coords_diffs(self.coords[frame])
        energies = self.energies[frame]
        self.energies_plot.set_xdata(coords_diffs)
        self.energies_plot.set_ydata(energies)

        return self.images, self.quiv, self.energies

    def animate(self):
        cycles = range(self.optimizer.cur_cycle)
        self.animation = animation.FuncAnimation(self.fig,
                                                 self.func,
                                                 frames=cycles,
                                                 interval=500)
        plt.show()

    def on_click(self, event):
        """Pause on any keypress."""
        #https://stackoverflow.com/questions/41557578
        if self.pause:
            self.animation.event_source.stop()
        else:
            self.animation.event_source.start()
        self.pause = not self.pause
