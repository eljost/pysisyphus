#!/usr/bin/env python3

from collections import namedtuple

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry

import matplotlib.pyplot as plt
import numpy as np


def get_geoms(coords=None):
    # if coords is None:
        # initial = np.array((-1.05274, 1.02776, 0))
        # final = np.array((1.94101, 3.85427, 0))
        # coords = (initial, final)
    if coords is None:
        left = np.array((0.188646, 1.45698, 0))
        right = np.array((0.950829, 1.54153, 0))
        left = np.array((0.354902, 1.34229, 0))
        right = np.array((0.881002, 1.71074, 0))
        # left = np.array((0.531642, 1.41899, 0))
        # right = np.array((0.702108, 1.57077, 0))
        # coords = (left, right)
        coords = (right, left)
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
    dimer_method(geoms, calc_getter)


DimerCycle = namedtuple("DimerCycle",
                        "org_coords trial_coords rot_coords trans_coords f0 f0_mod",
)


def perpendicular_force(force, vec):
    return force - force.dot(vec)*vec


def rotate_R1(coords, angle, N, theta, dR):
    """Only valid for rotation of R1!"""
    return coords + (N*np.cos(angle) + theta*np.sin(angle)) * dR


def dimer_method(geoms, calc_getter):
    """1 oder 2 startpunkte.
    bei 2 startpunkten dR bestimmen, be 1 startpunkt muss N zufällig? bestimmt
    werden und dR gegeben sein.

    Params
        trial_angle
        tiral_angle_base
        rot_trials              Max. anzahl rotationsversuche
        rot_force_thresh
    """

    geom1, geom2 = geoms 

    # N points from geom2 to geom1
    N = geom1.coords - geom2.coords
    N /= np.linalg.norm(N)

    dR = np.linalg.norm(geom1.coords - geom2.coords) / 2

    geom0 = geom2.copy()
    geom0.coords += dR*N
    geom0.set_calculator(calc_getter())

    # Parameters
    angle_base = np.deg2rad(5)
    step_base = 0.05

    dimer_cycles = list()
    for i in range(25):
        coords1 = geom1.coords
        coords0 = geom0.coords
        coords2 = geom2.coords

        # Forces at geom1 and geom2
        f1 = geom1.forces
        f2 = geom2.forces

        print(f"Cycle {i}")
        C = (f2 - f1).dot(N) / (2*dR)
        print(f"Curvature: {C:.4f}")

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.
        f1_perp = perpendicular_force(f1, N)
        f2_perp = perpendicular_force(f2, N)
        f_perp = f1_perp - f2_perp
        theta = f_perp / np.linalg.norm(f_perp)
        for i in range(5):
            angle = angle_base * (i+1)
            # print(f"{i}: {angle}")
            coords1_star = rotate_R1(coords0, angle, N, theta, dR)
            N_star = coords1_star - coords0
            N_star /= np.linalg.norm(N_star)
            coords2_star = coords0 - N_star*dR
            rg1 = geom0.copy()
            rg1.set_calculator(calc_getter())
            rg1.coords = coords1_star
            rg2 = geom0.copy()
            rg2.set_calculator(calc_getter())
            rg2.coords = coords2_star
            fr1 = rg1.forces
            fr2 = rg2.forces
            nr1 = np.linalg.norm(fr1)
            nr2 = np.linalg.norm(fr2)
            trial_coords = np.array((coords1_star, coords0, coords2_star))
            # print(f"\t forces1: {fr1}, {nr1}")
            # print(f"\t forces2: {fr1}, {nr2}")
            # print("Breaking after first rotation!")
            # print()
            break
        rot_diff = (fr1 - fr2).dot(theta)
        org_diff = (f1 - f2).dot(theta)
        F = (rot_diff + org_diff) / 2
        F_dash = (rot_diff - org_diff) / angle
        # angle_min = -0.5*np.arctan2(2*F, F_dash) - angle/2
        angle_min = -0.5*np.arctan2(F_dash, 2*F) - angle/2
        # print("angle", angle, "angle_min", angle_min)
        # print("angle", np.rad2deg(angle), "angle_min", np.rad2deg(angle_min))
        coords1_rot = rotate_R1(coords0, angle_min, N, theta, dR)
        hess = geom0.hessian
        w, v = np.linalg.eig(hess)
        trans_mode = v[:,0]
        # coords1_rot = rotate_R1(coords0, angle_base, N, theta, dR)
        N_rot = coords1_rot - coords0
        N_rot /= np.linalg.norm(N_rot)
        mode_ovlp = trans_mode.dot(N_rot)
        coords2_rot = coords0 - N_rot*dR
        geom1.coords = coords1_rot
        geom2.coords = coords2_rot

        # Translation
        # f0 = (f1 + f2) / 2
        f0 = geom0.forces
        nr0 = np.linalg.norm(f0)
        if C > 0:
            f0_mod = -f0.dot(N_rot) * N_rot
        else:
            f0_mod = f0 - 2*f0.dot(N_rot) * N_rot
        N0_mod = f0_mod / np.linalg.norm(f0_mod)

        # Small displacement of midpoint
        geom0_mod = geom0.copy()
        trans_midpoint_coords = coords0 + step_base * f0_mod


        coords1_trans = trans_midpoint_coords + dR*N_rot
        coords2_trans = trans_midpoint_coords - dR*N_rot


        org_coords = np.array((coords1, coords0, coords2))
        rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        trans_coords = np.array((coords1_trans, trans_midpoint_coords, coords2_trans))

        dc = DimerCycle(org_coords, trial_coords, rot_coords, trans_coords, f0, f0_mod)
        dimer_cycles.append(dc)

        geom1.coords = coords1_trans
        geom2.coords = coords2_trans
        geom0.coords = trans_midpoint_coords
        N = N_rot
        print(f"Overlap: {mode_ovlp:.4f}")
        print(f"norm(f0): {nr0:.4f}")
        print()

    # plot_geoms(dimer_cycles)
    plot_dimer_cycles(dimer_cycles)


if __name__ == "__main__":
    run()
