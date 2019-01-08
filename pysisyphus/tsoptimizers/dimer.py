#!/usr/bin/env python3

from collections import namedtuple

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry

import matplotlib.pyplot as plt
import numpy as np


def get_geoms(coords=None):
    if coords is None:
        left = np.array((0.188646, 1.45698, 0))
        right = np.array((0.950829, 1.54153, 0))
        left = np.array((0.354902, 1.34229, 0))
        right = np.array((0.881002, 1.71074, 0))
        left = np.array((0.531642, 1.41899, 0))
        right = np.array((0.702108, 1.57077, 0))
        coords = (right, left)
        # near_ts = np.array((0.553726, 1.45458, 0))
        # coords = (near_ts, )
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    for geom in geoms:
        geom.set_calculator(AnaPot())
    return geoms


def plot_dimer(dimer, ax, label=None, color=None, marker="o"):
    lines = ax.plot(dimer[:,0], dimer[:,1], marker=marker,
                    label=label, color=color)
    return lines


def plot_dimer_cycles(dimer_cycles):
    pot = AnaPot()
    levels = np.linspace(-2.8, 3.6, 50)
    pot.plot(levels=levels)

    ax = pot.ax
    cycs = len(dimer_cycles)
    last = 5
    for i, dc in enumerate(dimer_cycles[-last:]):
        cyc = cycs - last + i
        label = f"Cycle {cyc}"
        org_lbl = f"Org {cyc}"
        trial_lbl = f"Trial {cyc}"
        rot_lbl = f"Rot {cyc}"
        # org = plot_dimer(dc.org_coords, ax, label=org_lbl)
        # color = org[0].get_color()
        # trial = plot_dimer(dc.trial_coords, ax,
                         # label=trial_lbl, color=color, marker="x")
        # rot = plot_dimer(dc.rot_coords, ax,
                         # label=rot_lbl, color=color, marker=".")
        rot = plot_dimer(dc.rot_coords, ax,
                         label=rot_lbl, marker=".")
    pot.ax.legend()
    plt.show()


def plot_geoms(dimer_cycles):
    coords = np.array([geom.coords for geom in geoms])
    pot = AnaPot()
    levels = np.linspace(-2.8, 3.6, 50)
    pot.plot(levels=levels)
    xs = coords[:,0]
    ys = coords[:,1]
    pot.ax.scatter(coords[:,0], coords[:,1])
    for i, (x, y) in enumerate(zip(xs, ys)):
        pot.ax.scatter(x, y, label=f"{i}")
    pot.ax.legend()
    plt.show()


def run():
    geoms = get_geoms()
    calc_getter = AnaPot
    dimer_cycles = dimer_method(geoms, calc_getter, ana_2dpot=True)
    plot_dimer_cycles(dimer_cycles)



DimerCycle = namedtuple("DimerCycle",
                        "org_coords trial_coords rot_coords trans_coords f0 f0_mod",
)


def make_unit_vec(vec1, vec2):
    """Return unit vector pointing from vec2 to vec1."""
    diff = vec1 - vec2
    return diff / np.linalg.norm(diff)


def perpendicular_force(force, vec):
    return force - force.dot(vec)*vec


def rotate_R1(coords, angle, N, theta, dR):
    """Only valid for rotation of R1!"""
    return coords + (N*np.cos(angle) + theta*np.sin(angle)) * dR


def get_curvature(f1, f2, N, dR):
    return (f2 - f1).dot(N) / (2*dR)


def get_geom_getter(ref_geom, calc_setter):
    def geom_from_coords(coords):
        new_geom = ref_geom.copy()
        new_geom.coords = coords
        new_geom.set_calculator(calc_setter())
        return new_geom
    return geom_from_coords


def dimer_method(geoms, calc_getter, ana_2dpot=False):
    """Dimer method using steepest descent for rotation and translation.

    See
        # Original paper
        [1] https://doi.org/10.1063/1.480097
        # Improved dimer
        [2] https://doi.org/10.1063/1.2104507
        # Several trial rotations
        [3] https://doi.org/10.1063/1.1809574
        # Superlinear dimer
        [4] https://doi.org/10.1063/1.2815812
    """
    # Parameters
    trial_angle = np.deg2rad(5)
    angle_thresh = np.deg2rad(5)
    trial_step = 0.05
    dR_base = 0.1
    max_cycles = 15

    geom_getter = get_geom_getter(geoms[0], calc_getter)

    assert len(geoms) in (1, 2), "geoms argument must be of length 1 or 2!"
    if len(geoms) == 2:
        geom1, geom2 = geoms
        dR = np.linalg.norm(geom1.coords - geom2.coords) / 2
        # N points from geom2 to geom1
        N = make_unit_vec(geom1.coords, geom2.coords)
        coords0 = geom2.coords + dR*N
        geom0 = geom_getter(coords0)
    # This handles cases where only one geometry is supplied. We use the
    # geometry as midpoint, select a random dimer direction and derive
    # geom1 and geom2 from it.
    else:
        geom0 = geoms[0]
        # Assign random unit vector and use dR_base to for dR
        coord_size = geom0.coords.size
        N = np.random.rand(coord_size)
        if ana_2dpot:
            N[2] = 0
        N /= np.linalg.norm(N)
        coords1 = geom0.coords + dR_base*N
        geom1 = geom_getter(coords1)
        coords2 = geom0.coords - dR_base*N
        geom2 = geom_getter(coords2)
        dR = dR_base

    dimer_cycles = list()
    for i in range(max_cycles):
        coords0 = geom0.coords
        coords1 = geom1.coords
        coords2 = geom2.coords

        # Forces at geom1 and geom2
        f1 = geom1.forces
        f2 = geom2.forces

        C = get_curvature(f1, f2, N, dR)

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.
        f1_perp = perpendicular_force(f1, N)
        f2_perp = perpendicular_force(f2, N)
        f_perp = f1_perp - f2_perp
        theta = f_perp / np.linalg.norm(f_perp)
        # Trial rotation for finite difference calculation of rotational force
        # and rotational curvature.
        coords1_star = rotate_R1(coords0, trial_angle, N, theta, dR)
        N_star = make_unit_vec(coords1_star, coords0)
        coords2_star = coords0 - N_star*dR
        geom1_star = geom_getter(coords1_star)
        geom2_star = geom_getter(coords2_star)
        f1_star = geom1_star.forces
        f2_star = geom2_star.forces
        trial_coords = np.array((coords1_star, coords0, coords2_star))

        rot_diff = (f1_star - f2_star).dot(theta)
        org_diff = (f1 - f2).dot(theta)
        # Rotational force
        F = (rot_diff + org_diff) / 2
        # Rotational curvature
        F_dash = (rot_diff - org_diff) / trial_angle
        angle_min = -0.5*np.arctan2(F_dash, 2*F) - trial_angle/2
        if angle_min > angle_thresh:
            coords1_rot = rotate_R1(coords0, angle_min, N, theta, dR)
            N_rot = make_unit_vec(coords1_rot, coords0)
            coords2_rot = coords0 - N_rot*dR
            geom1.coords = coords1_rot
            geom2.coords = coords2_rot
        else:
            # Don't do rotation for small angles
            N_rot = N
            coords1_rot = coords1
            coords2_rot = coords2
            print("angle_min below threshold! no rotation!")

        # Translation
        f0 = geom0.forces
        nr0 = np.linalg.norm(f0)
        f0_mod = -f0.dot(N_rot)*N_rot
        if C < 0:
            f0_mod = f0 + 2*f0_mod
        N0_mod = f0_mod / np.linalg.norm(f0_mod)

        # Small displacement of midpoint
        geom0_mod = geom0.copy()
        trans_midpoint_coords = coords0 + trial_step * f0_mod

        coords1_trans = trans_midpoint_coords + dR*N_rot
        coords2_trans = trans_midpoint_coords - dR*N_rot

        # Save cycle information
        org_coords = np.array((coords1, coords0, coords2))
        rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        trans_coords = np.array((coords1_trans, trans_midpoint_coords, coords2_trans))
        dc = DimerCycle(org_coords, trial_coords, rot_coords, trans_coords, f0, f0_mod)
        dimer_cycles.append(dc)

        # Update dimer coordinates for next cycle
        geom1.coords = coords1_trans
        geom2.coords = coords2_trans
        geom0.coords = trans_midpoint_coords
        N = N_rot

        print(f"Cycle {i}")
        print(f"Initial curvature: {C:.4f}")
        print(f"norm(f0): {nr0:.4f}")
        if ana_2dpot:
            hess = geom0.hessian
            w, v = np.linalg.eig(hess)
            trans_mode = v[:,0]
            mode_ovlp = trans_mode.dot(N_rot)
            print(f"Overlap: {mode_ovlp:.4f}")
        print()
    return dimer_cycles


if __name__ == "__main__":
    run()
