#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from calculators.MullerBrownPot import MullerBrownPot
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent
from optimizers.NaiveSteepestDescent import NaiveSteepestDescent

CYCLES = 2
IMAGES = 20

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


def get_geoms():
    educt = np.array((0.6215, 0.02838, 0)) # Minimum B
    product = np.array((-0.558, 1.442, 0)) # Minimum A
    #product = np.array((-0.05, 0.467, 0)) # Minimum C
    #product = np.array((-0.822, 0.624, 0)) # Saddle point
    atoms = ("H", )
    geoms = [Geometry(atoms, coords) for coords in (educt, product)]
    return geoms

def run_cos_opt(cos_class, reparametrize=False):
    geoms = get_geoms()
    cos = cos_class(geoms)
    cos.interpolate_images(IMAGES)
    for img in cos.images[1:-1]:
        img.set_calculator(MullerBrownPot())

    sd = NaiveSteepestDescent(cos,
                         max_cycles=CYCLES,
                         max_force_thresh=0.05,
                         rms_force_thresh=0.01,
                         alpha=-0.005)
    if reparametrize:
        sd.run(reparam=cos.reparametrize)
    else:
        sd.run()
    plot_cos_opt(sd)


if __name__ == "__main__":
    #run_cos_opt(NEB)
    run_cos_opt(SimpleZTS, reparametrize=True)
    #plot_mullerbrownpot()
    #plt.show()

