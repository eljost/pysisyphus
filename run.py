#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from cos.NEB import NEB
from cos.NoSpringNeb import NoSpringNEB
from cos.SimpleZTS import SimpleZTS

from calculators.AnaPot import AnaPot

from qchelper.geometry import parse_xyz_file
from Geometry import Geometry
from calculators.ORCA import ORCA
from optimizers.SteepestDescent import SteepestDescent


def plot_neb(neb):
    x = np.linspace(-2, 2.5, 100)
    y = np.linspace(0, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = neb.calculator.get_energy(X, Y)

    width = 8
    height = width

    fig, ax = plt.subplots(figsize=(width, height))

    # Potential
    levels = np.linspace(-4, 8, 20)
    contours = ax.contour(X, Y, Z, levels)
    ax.clabel(contours, inline=1, fontsize=10)

    img_xs = neb.old_images[:,0]
    img_ys = neb.old_images[:,1]

    # Image positions
    ax.plot(img_xs, img_ys, "ro", ls="-")

    # Gradient
    ax.quiver(img_xs, img_ys, neb.grad_xs, neb.grad_ys)

    # True force component
    #C = np.hypot(true_forces[:,0], true_forces[:,1])
    #ax.quiver(img_xs[1:-1], img_ys[1:-1], true_forces[:,0], true_forces[:,1], C)

    # Total force component
    C = np.hypot(neb.total_forces[:,0], neb.total_forces[:,1])
    ax.quiver(img_xs[1:-1], img_ys[1:-1], neb.total_forces[:,0], neb.total_forces[:,1], C)

    # New images
    ax.plot(neb.images[:,0], neb.images[:,1], "bx", ls="-")

    return fig


def plot_anapot_neb(calculator, optimizer):
    x = np.linspace(-2, 2.5, 100)
    y = np.linspace(0, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, 0)
    fake_atoms = ("H", "H")
    coords = np.stack((X, Y, Z))
    pot = calculator.get_energy(fake_atoms, coords)

    width = 8
    height = width

    fig, ax = plt.subplots(figsize=(width, height))

    # Potential
    levels = np.linspace(-4, 8, 20)
    contours = ax.contour(X, Y, pot, levels)
    ax.clabel(contours, inline=1, fontsize=10)
    plt.show()
    import sys; sys.exit()

    img_xs = neb.old_images[:,0]
    img_ys = neb.old_images[:,1]

    # Image positions
    ax.plot(img_xs, img_ys, "ro", ls="-")

    # Gradient
    ax.quiver(img_xs, img_ys, neb.grad_xs, neb.grad_ys)

    # True force component
    #C = np.hypot(true_forces[:,0], true_forces[:,1])
    #ax.quiver(img_xs[1:-1], img_ys[1:-1], true_forces[:,0], true_forces[:,1], C)

    # Total force component
    C = np.hypot(neb.total_forces[:,0], neb.total_forces[:,1])
    ax.quiver(img_xs[1:-1], img_ys[1:-1], neb.total_forces[:,0], neb.total_forces[:,1], C)

    # New images
    ax.plot(neb.images[:,0], neb.images[:,1], "bx", ls="-")

    return fig

"""
def run_new():
    neb.save("iter000.trj")
    ics = list()
    for i in range(1, 7):
        neb.cycle()
        inner_coords = neb.images[1].coords
        ics.append(inner_coords)
        neb.save("iter{:03}.trj".format(i))

    ic_reshaped = [ic.reshape((-1,3)) for ic in ics]
    atoms = atoms_coords[0][0]
    from qchelper.geometry import make_trj_str
    trj_str = make_trj_str(atoms, ic_reshaped)
    with open("inner_coords.trj", "w") as handle:
        handle.write(trj_str)
"""

def run_neb():
    educt = "h2o_inv_educt.xyz"
    product = "h2o_inv_product.xyz" 
    xyz_fns = [educt, product]
    atoms_coords = [parse_xyz_file(fn) for fn in xyz_fns]
    geoms = [Geometry(atoms, coords.flatten()) for atoms, coords in atoms_coords]
    neb = NEB(geoms)
    neb.interpolate_images(2)
    for img in neb.images[1:-1]:
        img.set_calculator(ORCA())

    steepest_descent(neb)


def run_anapot_neb():
    educt = np.array((-1.05274, 1.02776, 0))
    product = np.array((1.94101, 3.85427, 0))
    atoms = ("H", "H")
    geoms = [Geometry(atoms, coords) for coords in (educt, product)]
    neb = NEB(geoms)
    neb.interpolate_images(7)
    for img in neb.images[1:-1]:
        img.set_calculator(AnaPot())

    sd = SteepestDescent(neb)#, alpha=-0.05, max_cycles=20)
    sd.cycle()
    fig = plot_anapot_neb(AnaPot(), sd)


def run_opt():
    atoms, coords = parse_xyz_file("h2o_inv_educt.xyz")
    geoms = Geometry(atoms, coords.flatten())
    geoms.set_calculator(ORCA())
    steepest_descent(geoms)

if __name__ == "__main__":
    #run_neb()
    #run_opt()
    run_anapot_neb()

