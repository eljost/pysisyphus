#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from calculators.AnaPot import AnaPot
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent

CYCLES = 5
IMAGES = 7

def plot_anapot():
    x = np.linspace(-2, 2.5, 100)
    y = np.linspace(0, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, 0)
    fake_atoms = ("H", "H")
    pot_coords = np.stack((X, Y, Z))
    pot = AnaPot().get_energy(fake_atoms, pot_coords)["energy"]

    width = 8
    height = width

    fig, ax = plt.subplots(figsize=(width, height))

    # Potential
    levels = np.linspace(-4, 8, 20)
    contours = ax.contour(X, Y, pot, levels)
    ax.clabel(contours, inline=1, fontsize=10)

    return fig, ax


def plot_anapot_opt(optimizer):
    coords = optimizer.coords
    forces = optimizer.forces
    coords = [c.reshape((-1, 3)) for c in coords]
    forces = [f.reshape((-1, 3)) for f in forces]

    for cycle in range(optimizer.cur_cycle):
        fig, ax = plot_anapot()
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
    educt = np.array((-1.05274, 1.02776, 0))
    product = np.array((1.94101, 3.85427, 0))
    atoms = ("H", "H")
    geoms = [Geometry(atoms, coords) for coords in (educt, product)]
    return geoms


def run_anapot_neb():
    geoms = get_geoms()
    neb = NEB(geoms)
    neb.interpolate_images(IMAGES)
    for img in neb.images[1:-1]:
        img.set_calculator(AnaPot())

    sd = SteepestDescent(neb,
                         max_cycles=CYCLES,
                         max_force_thresh=0.05,
                         rms_force_thresh=0.01,
                         alpha=-0.05)
    sd.run()
    plot_anapot_opt(sd)


def run_anapot_szts():
    geoms = get_geoms()
    szts = SimpleZTS(geoms)
    szts.interpolate_images(IMAGES)
    for img in szts.images[1:-1]:
        img.set_calculator(AnaPot())
    sd = SteepestDescent(szts,
                         max_cycles=CYCLES,
                         max_force_thresh=0.05,
                         rms_force_thresh=0.01,
                         alpha=-0.05)
    #sd.run()
    sd.run(reparam=szts.reparametrize)
    #plot_anapot_opt(sd)


if __name__ == "__main__":
    #run_anapot_neb()
    run_anapot_szts()

