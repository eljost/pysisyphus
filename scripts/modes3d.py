#!/usr/bin/env python3

import argparse
import itertools as it
from math import log
import sys

import numpy as np
from ursina import (
    color,
    Cylinder,
    destroy,
    DirectionalLight,
    EditorCamera,
    Entity,
    Func,
    Sequence,
    Text,
    Ursina,
    window,
)
from ursina.shaders import lit_with_shadows_shader

from pysisyphus.helpers import geom_loader
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


def get_bond_data(coords3d, bonds):
    directions = list()
    vals = list()
    xyz = list()
    for bond in bonds:
        from_, to_ = bond
        direction = coords3d[to_] - coords3d[from_]
        direction /= np.linalg.norm(direction)
        direction = tuple(direction)
        val = Stretch._calculate(coords3d, bond)
        x, y, z = coords3d[from_]
        directions.append(direction)
        vals.append(val)
        xyz.append((x, y, z))
    return directions, vals, xyz


def get_bond_data_from_geom(geom):
    return get_bond_data(geom.coords3d, geom.internal.bond_atom_indices)


def render_atoms(atoms, coords3d):
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
    return spheres


def render_bonds(coords3d, bonds=None):
    if bonds is None:
        bonds = list()

    directions, vals, xyz = get_bond_data(coords3d, bonds)
    cylinders = list()
    for direction, val, (x, y, z) in zip(directions, vals, xyz):
        cyl_ = Cylinder(direction=direction, radius=0.1, height=val)
        cyl = Entity(
            model=cyl_,
            x=x,
            y=y,
            z=z,
            shader=lit_with_shadows_shader,
        )
        cyl.hide()
        cylinders.append(cyl)
    return cylinders


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


def palindrome_cycler(iterable):
    inds = list(range(len(iterable)))
    cycler = it.cycle(inds + inds[1:-1][::-1])
    return cycler


def from_h5_hessian(fn):
    geom = geom_from_hessian(fn, coord_type="redund")
    geoms = (geom,)

    # Cartesian displacements
    nus, *_, displs = geom.get_normal_modes()
    nus[1] = -2000
    displs[:, 1] = displs[:, -3]
    # Imaginary frequencies
    imag_mask = nus < 0
    imag_trj_coords = list()
    for i, nu in enumerate(nus[imag_mask]):
        print(f"Imaginary Mode {i:02d}: {nu:.4f} cm⁻¹")
        imag_trj_coords.append(get_tangent_trj_coords(geom.cart_coords, displs[:, i]))

    imag_num = imag_mask.sum()
    global imag_cycler
    global imag_ind
    imag_cycler = it.cycle(range(imag_num))
    imag_ind = next(imag_cycler)
    c3d_cycler = palindrome_cycler(imag_trj_coords[0])
    return geoms, c3d_cycler, imag_trj_coords


def from_geom(fn):
    geoms = geom_loader(fn, coord_type="redund", iterable=True)
    trj_coords = [geom.coords3d for geom in geoms]
    cycler = palindrome_cycler(trj_coords)
    return geoms, cycler, trj_coords


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fn")
    parser.add_argument("--delay", type=float, default=0.025)
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    fn = args.fn
    delay = args.delay

    if fn.endswith(".h5"):
        geoms, cycler, trj_coords = from_h5_hessian(fn)
    else:
        geoms, cycler, trj_coords = from_geom(fn)

    window.title = "pysisyphus molecule viewer"
    window.borderless = False
    window.windowed_resolution = (1920, 1080)

    app = Ursina()
    window.update_aspect_ratio()
    EditorCamera(rotate_around_mouse_hit=True)

    geom = geoms[0]
    spheres = render_atoms(geom.atoms, geom.coords3d)
    cylinders = [
        render_bonds(geom.coords3d, geom.internal.bond_atom_indices) for geom in geoms
    ]
    cur_cylinders = cylinders[0]

    def update_positions():
        nonlocal cur_cylinders

        ind = next(cycler)
        try:
            c3d = trj_coords[imag_ind][ind]
            # No bond update yet for Hessian
            new_cylinders = cur_cylinders
        except NameError:
            c3d = trj_coords[ind]
            new_cylinders = cylinders[ind]

        for sphere, xyz in zip(spheres, c3d):
            sphere.position = xyz

        for cylinder in cur_cylinders:
            cylinder.hide()

        for cylinder in new_cylinders:
            cylinder.show()

        cur_cylinders = new_cylinders

    # Light and shadows
    pivot = Entity()
    DirectionalLight(parent=pivot, shadows=True)
    seq = Sequence(delay, Func(update_positions), loop=True)
    seq.start()
    app.run()


if __name__ == "__main__":
    run()
