#!/usr/bin/env python3

import argparse
import itertools as it
from math import log
import sys

import numpy as np
from ursina import *
from ursina.shaders import lit_with_shadows_shader

from pysisyphus.intcoords.Stretch import Stretch
from pysisyphus.io.hessian import geom_from_hessian


C = {
    "h": color.white,
    "c": color.gray,
    "n": color.blue,
    "o": color.red,
    "f": color.lime,
    "s": color.yellow,
    "cl": color.green,
    "p": color.orange,
}
S = {
    "h": 1.0,
    "n": 1.25,
    "o": 1.25,
    "s": 1.25,
}


def render_molecule(atoms, coords3d, bonds=None):
    if bonds is None:
        bonds = list()

    spheres = list()
    for atom, (x, y, z) in zip(atoms, coords3d):
        atom = atom.lower()
        scale_ = S.get(atom, 1.0)
        scale = (scale_, scale_, scale_)
        acolor = C.get(atom, color.pink)
        sphere = Entity(
            model="sphere",
            color=acolor,
            scale=scale,
            world_x=x,
            world_y=y,
            world_z=z,
            shader=lit_with_shadows_shader,
        )
        spheres.append(sphere)

    bonds_ = list()
    for i, bond in enumerate(bonds):
        from_, to_ = bond
        direction = coords3d[to_] - coords3d[from_]
        direction /= np.linalg.norm(direction)
        direction = tuple(direction)
        val = Stretch._calculate(coords3d, bond)
        x, y, z = coords3d[from_]

        cyl_ = Cylinder(direction=direction, radius=0.1, height=val)
        cyl = Entity(
            model=cyl_,
            x=x,
            y=y,
            z=z,
            shader=lit_with_shadows_shader,
        )
        bonds_.append(cyl)
    return spheres, bonds


def input(key):
    global imag_ind
    global imag_cycler

    if key == "space":
        imag_ind = next(imag_cycler)
        msg = f"Root {imag_ind}"
        popup_text = Text(msg, x=0.30, y=0.45)
        destroy(popup_text, delay=0.5)


def get_tangent_trj_coords(coords, tangent, points=10, displ=None):
    atom_num = coords.reshape(-1, 3).shape[0]
    if displ is None:
        displ = 0.43429 * log(atom_num) + 0.52288
    step_sizes = np.linspace(-displ, displ, 2 * points + 1)
    steps = step_sizes[:, None] * tangent
    trj_coords = coords[None, :] + steps
    trj_coords = trj_coords.reshape(step_sizes.size, -1, 3)
    return trj_coords


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("hessian")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    geom = geom_from_hessian(args.hessian, coord_type="redund")

    # Cartesian displacements
    nus, *_, displs = geom.get_normal_modes()
    nus[1] = -2000
    displs[:, 1] = displs[:, -3]
    # Imaginary frequencies
    imag_mask = nus < 0
    imag_trj_coords = list()
    for i, (nu, displ) in enumerate(zip(nus[imag_mask], displs.T[imag_mask])):
        print(f"Imaginary Mode {i:02d}: {nu:.4f} cm⁻¹")
        imag_trj_coords.append(get_tangent_trj_coords(geom.cart_coords, displs[:, i]))

    imag_num = imag_mask.sum()
    global imag_cycler
    global imag_ind
    imag_cycler = it.cycle(range(imag_num))
    imag_ind = next(imag_cycler)
    inds = list(range(len(imag_trj_coords[0])))
    # Palindrome
    c3d_cycler = it.cycle(inds + inds[1:-1][::-1])

    def update_positions():
        ind = next(c3d_cycler)
        c3d = imag_trj_coords[imag_ind][ind]
        for sphere, xyz in zip(spheres, c3d):
            sphere.position = xyz

    window.title = "pysisyphus normal mode viewer"
    window.borderless = False
    window.windowed_resolution = (1920, 1080)
    app = Ursina()
    window.update_aspect_ratio()
    EditorCamera(rotate_around_mouse_hit=True)
    spheres, bonds = render_molecule(
        geom.atoms, geom.coords3d, geom.internal.bond_atom_indices
    )
    # Light and shadows
    pivot = Entity()
    DirectionalLight(parent=pivot, shadows=True)
    seq = Sequence(0.05, Func(update_positions), loop=True)
    seq.start()
    app.run()


if __name__ == "__main__":
    run()
