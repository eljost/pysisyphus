#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class AnimPlot:

    def __init__(self, calculator, optimizer):
        self.calculator = calculator
        self.optimizer = optimizer

        coords = self.optimizer.coords
        coords = [c.reshape(-1, 3) for c in coords]
        self.coords = coords

        self.fig, self.ax = plt.subplots()


    #def init_func(self):
        x = np.linspace(-1.75, 1.25, 100)
        y = np.linspace(-0.5, 2.25, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y, Z))
        pot = self.calculator().get_energy(fake_atoms, pot_coords)["energy"]

        levels = np.linspace(-150, 5, 30)
        contours = self.ax.contour(X, Y, pot, levels)
        self.ax.clabel(contours, inline=1, fontsize=5)

        imagex = self.coords[0][:,0]
        imagey = self.coords[0][:,1]
        self.images, = self.ax.plot(imagex, imagey, "ro", ls="-")
        #p = self.ax.plot(imagex, imagey, "ro", ls="-")

    def func(self, frame):
        self.fig.suptitle("Cycle {}".format(frame))
        imagex = self.coords[frame][:,0]
        imagey = self.coords[frame][:,1]
        self.images.set_xdata(imagex)
        self.images.set_ydata(imagey)
        #= self.ax.plot(imagex, imagey, "ro", ls="-")
        return self.images, 

    def animate(self):
        cycles = range(1, self.optimizer.cur_cycle)
        an = animation.FuncAnimation(self.fig,
                                     self.func,
                                     frames=cycles,
                                     #init_func=self.init_func,
                                     interval=750)
        plt.show()

    """
    def plot_mullerbrownpot():
        x = np.linspace(-1.75, 1.25, 100)
        y = np.linspace(-0.5, 2.25, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y, Z))
        pot = MullerBrownPot().get_energy(fake_atoms, pot_coords)["energy"]

        width = 8
        height = width

        fig, ax = plt.subplots(figsize=(width, height))

        # Potential
        levels = np.linspace(-150, 5, 30)
        contours = ax.contour(X, Y, pot, levels)
        ax.clabel(contours, inline=1, fontsize=5)

        return fig, ax


    def plot_cos_opt(optimizer):
        coords = optimizer.coords
        forces = optimizer.forces
        coords = [c.reshape((-1, 3)) for c in coords]
        forces = [f.reshape((-1, 3)) for f in forces]

        for cycle in range(optimizer.cur_cycle):
            fig, ax = plot_mullerbrownpot()
            fig.suptitle("Cycle {}".format(cycle))

            imagex = coords[cycle][:,0]
            imagey = coords[cycle][:,1]
            ax.plot(imagex, imagey, "ro", ls="-")

            # Force
            forcesx = forces[cycle][:,0]
            forcesy = forces[cycle][:,1]
            ax.quiver(imagex, imagey, forcesx, forcesy)

            #plt.tight_layout()
            plt.show()
    """
