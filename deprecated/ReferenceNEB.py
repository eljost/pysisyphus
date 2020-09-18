#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from sympy import symbols, diff, lambdify, sympify


class AnimPlot:

    def __init__(self, coords, tangents, cycles, ge):
        self.coords = coords
        self.tangents = tangents
        print("tangents hsape", tangents.shape)
        self.cycles = cycles
        self.ge = ge

        self.fig, self.ax = plt.subplots(figsize=(8,8))
        self.pause = True
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

        xlim = (-2, 2.5)
        ylim = (0, 5)
        x = np.linspace(*xlim, 100)
        y = np.linspace(*ylim, 100)
        X, Y = np.meshgrid(x, y)
        pot_coords = np.stack((X, Y))
        pot = self.ge(pot_coords)

        levels = (-3, 6, 50)
        levels = np.linspace(*levels)
        contours = self.ax.contour(X, Y, pot, levels)

        images_x = self.coords[0][:,0]
        images_y = self.coords[0][:,1]
        """
        forces_x = self.forces[0][:,0]
        forces_y = self.forces[0][:,1]
        """
        tangents_x = self.tangents[0][:,0]
        tangents_y = self.tangents[0][:,1]

        self.images, = self.ax.plot(images_x, images_y, "ro", ls="-")
        # Tangents
        self.tangent_quiv = self.ax.quiver(images_x[1:-1], images_y[1:-1],
                                           tangents_x, tangents_y, color="b")


    def func(self, frame):
        self.fig.suptitle("Cycle {}".format(frame))

        images_x = self.coords[frame][:,0]
        images_y = self.coords[frame][:,1]
        self.images.set_xdata(images_x)
        self.images.set_ydata(images_y)

        offsets = np.stack((images_x[1:-1], images_y[1:-1]), axis=-1).flatten()

        # Update tangent quiver
        tangents_x = self.tangents[frame][:,0]
        tangents_y = self.tangents[frame][:,1]
        self.tangent_quiv.set_offsets(offsets)
        self.tangent_quiv.set_UVC(tangents_x, tangents_y)

    def animate(self):
        self.animation = animation.FuncAnimation(self.fig,
                                                 self.func,
                                                 frames=range(self.cycles),
                                                 interval=250)
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


def make_potential(V_str):
    x, y = symbols("x y")
    V = sympify(V_str)
    dVdx = diff(V, x)
    dVdy = diff(V, y)
    V = lambdify((x, y), V, "numpy")
    dVdx = lambdify((x, y), dVdx, "numpy")
    dVdy = lambdify((x, y), dVdy, "numpy")
    return V, dVdx, dVdy


def get_energy(V, coords):
    x, y = coords
    #x = coords[:, 0]
    #y = coords[:, 1]
    return V(x, y)


def get_forces(dVdx, dVdy, coords):
    x, y = coords
    dVdx = dVdx(x, y)
    dVdy = dVdy(x, y)
    return -np.array((dVdx, dVdy))


def interpolate(initial, final, new_images):
    step = (final-initial) / (new_images+1)
    # initial + i*step
    i_array = np.arange(0, new_images+2)
    new_coords = initial + i_array[:, None]*step
    return np.array(new_coords)


def run():
    new_images = 8
    max_cycles = 10
    k = 0.01
    alpha = 0.05

    #V_str = "(1 - x**2 - y**2)**2 + (y**2) / (x**2 + y**2)"
    #initial = np.array((-0.5, 0.5, 0))
    #final = np.array((0.5, 0.5, 0))

    V_str = "4 + 4.5*x - 4*y + x**2 + 2*y**2-2*x*y + x**4 - 2*x**2*y"
    initial = np.array((-1.05274, 1.02776))
    final = np.array((1.94101, 3.85427))

    V, dVdx, dVdy = make_potential(V_str)
    ge = lambda c: get_energy(V, c)
    gf = lambda c: get_forces(dVdx, dVdy, c)

    coords = interpolate(initial, final, new_images=new_images)

    all_coords = list()
    all_true_forces = list()
    all_neb_forces = list()
    all_energies = list()
    all_tangents = list()
    for i in range(max_cycles):
        energies = ge(coords.transpose())
        forces = gf(coords.transpose())

        all_coords.append(coords.copy())
        all_true_forces.append(forces)
        all_energies.append(energies)

        neb_forces = list()
        tangents = list()
        for j in range(1, len(coords)-1):
            prev_coords = coords[j-1]
            jth_coords = coords[j]
            next_coords = coords[j+1]

            prev_energy = energies[j-1]
            jth_energy = energies[j]
            next_energy = energies[j+1]
            if next_energy > jth_energy > prev_energy:
                tangent = next_coords - jth_coords
            elif next_energy < jth_energy < prev_energy:
                tangent = jth_coords - prev_coords
            else:
                max_energy_diff = max(abs(next_energy-jth_energy),
                                      abs(prev_energy-jth_energy)
                )
                min_energy_diff = min(abs(next_energy-jth_energy),
                                      abs(prev_energy-jth_energy)
                )
                left_coords_diff = next_coords - jth_coords
                right_coords_diff = jth_coords - prev_coords

                if next_energy > prev_energy:
                    tangent = (left_coords_diff*max_energy_diff
                               + right_coords_diff*min_energy_diff
                    )
                else:
                    tangent = (left_coords_diff*min_energy_diff
                               + right_coords_diff*max_energy_diff
                    )
            tangent = tangent / np.linalg.norm(tangent)
            tangents.append(tangent)

            parallel_forces = (k * (np.linalg.norm(next_coords-jth_coords)
                                  - np.linalg.norm(prev_coords-jth_coords))
                               * tangent
            )
            perpendicular_gradient = (-forces[:,j]
                                      + np.dot(forces[:,j], tangent) * tangent
            )
            neb_forces.append(parallel_forces-perpendicular_gradient)
        neb_forces = np.array(neb_forces)
        all_neb_forces.append(neb_forces)
        all_tangents.append(tangents)

        steps = alpha*neb_forces
        coords[1:-1] += steps

    all_coords = np.array(all_coords)
    all_tangents = np.array(all_tangents)

    anim_plot = AnimPlot(all_coords, all_tangents, max_cycles, ge)
    anim_plot.animate()


if __name__ == "__main__":
    run()
