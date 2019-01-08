#!/usr/bin/env python3

from collections import namedtuple

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.TablePrinter import TablePrinter

import matplotlib.pyplot as plt
import numpy as np


def get_geoms(coords=None):
    if coords is None:
        # left = np.array((0.188646, 1.45698, 0))
        # right = np.array((0.950829, 1.54153, 0))
        # left = np.array((0.354902, 1.34229, 0))
        # right = np.array((0.881002, 1.71074, 0))
        left = np.array((0.531642, 1.41899, 0))
        right = np.array((0.702108, 1.57077, 0))
        coords = (right, left)

        # near_ts = np.array((0.553726, 1.45458, 0))
        # coords = (near_ts, )

        # left_far = np.array((-0.455116, 0.926978, 0))
        # right_far = np.array((-0.185653, 1.02486, 0))
        # coords = (left_far, right_far)

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
    for i, dc in enumerate(dimer_cycles):
        label = f"Cycle {i}"
        org_lbl = f"Org {i}"
        trial_lbl = f"Trial {i}"
        rot_lbl = f"Rot {i}"
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


def run():
    geoms = get_geoms()
    calc_getter = AnaPot
    dimer_kwargs = {
        "max_step": 0.1,
        "ana_2dpot": True,
    }
    dimer_cycles = dimer_method(geoms, calc_getter, **dimer_kwargs)
    plot_dimer_cycles(dimer_cycles[-5:])


DimerCycle = namedtuple("DimerCycle",
                        "org_coords trial_coords rot_coords trans_coords f0 f0_mod",
)


def make_unit_vec(vec1, vec2):
    """Return unit vector pointing from vec2 to vec1."""
    diff = vec1 - vec2
    return diff / np.linalg.norm(diff)


def perpendicular_force(force, vec):
    return force - force.dot(vec)*vec


def rotate_R1(coords, rad, N, theta, dR):
    """Only valid for rotation of R1!"""
    return coords + (N*np.cos(rad) + theta*np.sin(rad)) * dR


def get_curvature(f1, f2, N, dR):
    return (f2 - f1).dot(N) / (2*dR)


def get_f_mod(f, N, C):
    f_mod = -f.dot(N)*N
    if C < 0:
        f_mod = f + 2*f_mod
    return f_mod


def get_geom_getter(ref_geom, calc_setter):
    def geom_from_coords(coords):
        new_geom = ref_geom.copy()
        new_geom.coords = coords
        new_geom.set_calculator(calc_setter())
        return new_geom
    return geom_from_coords


def dimer_method(geoms, calc_getter,
                 max_step=0.3, max_cycles=50,
                 trial_angle=5, angle_thresh=5, ana_2dpot=False):
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

        To add:
            Comparing curvatures and adding π/2 if appropriate.
    """
    # Parameters
    max_cycles = 50
    dR_base = 0.1
    trial_rad = np.deg2rad(trial_angle)
    angle_thresh_rad = np.deg2rad(angle_thresh)
    dx = 0.05

    header = "Cycle Curvature max(f0) rms(f0)".split()
    col_fmts = "int float float float".split()
    table = TablePrinter(header, col_fmts)

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

    table.print_header()
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
        coords1_star = rotate_R1(coords0, trial_rad, N, theta, dR)
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
        F_dash = (rot_diff - org_diff) / trial_rad
        rad_min = -0.5*np.arctan2(F_dash, 2*F) - trial_rad/2
        if rad_min > angle_thresh_rad:
            coords1_rot = rotate_R1(coords0, rad_min, N, theta, dR)
            N_rot = make_unit_vec(coords1_rot, coords0)
            coords2_rot = coords0 - N_rot*dR
            geom1.coords = coords1_rot
            geom2.coords = coords2_rot
        else:
            # Don't do rotation for small angles
            N_rot = N
            coords1_rot = coords1
            coords2_rot = coords2
            table.print("Rotation angle too small. Skipping rotation.")

        # Translation
        f0 = geom0.forces
        f0_rms = np.sqrt(np.power(f0, 2).mean())
        f0_max = f0.max()

        row_args = [i, C, f0_rms, f0_max]
        table.print_row(row_args)

        if f0_rms <= 1e-3 and f0_max <= 1.5e-3:
            table.print("Converged!")
            break
        f0_mod = -f0.dot(N_rot)*N_rot
        f0_mod = get_f_mod(f0, N_rot, C)
        N0_mod = f0_mod / np.linalg.norm(f0_mod)
        trial_step = N0_mod * dx

        trial_coords0 = coords0 + trial_step*f0_mod
        trial_geom0 = geom_getter(trial_coords0)
        trial_f0 = trial_geom0.forces
        trial_f0_mod = get_f_mod(trial_f0, N_rot, C)
        F_mod = (trial_f0_mod + f0_mod).dot(N0_mod) / 2
        C_mod = (trial_f0_mod - f0_mod).dot(N0_mod) / dx
        step = (-F_mod/C_mod + dx/2)*N0_mod
        if np.linalg.norm(step) > max_step:
            step = max_step * f0_mod

        # Small displacement of midpoint
        coords0_trans = coords0 + step
        coords1_trans = coords0_trans + dR*N_rot
        coords2_trans = coords0_trans - dR*N_rot

        # Save cycle information
        org_coords = np.array((coords1, coords0, coords2))
        rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        trans_coords = np.array((coords1_trans, coords0_trans, coords2_trans))
        dc = DimerCycle(org_coords, trial_coords, rot_coords, trans_coords, f0, f0_mod)
        dimer_cycles.append(dc)

        # Update dimer coordinates for next cycle
        geom1.coords = coords1_trans
        geom2.coords = coords2_trans
        geom0.coords = coords0_trans
        N = N_rot

        header = "Cycle Curvature max(f0) rms(f0)".split()
        if ana_2dpot:
            hess = geom0.hessian
            w, v = np.linalg.eig(hess)
            trans_mode = v[:,0]
            mode_ovlp = trans_mode.dot(N_rot)
    return dimer_cycles


if __name__ == "__main__":
    run()
