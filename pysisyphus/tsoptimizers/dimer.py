#!/usr/bin/env python3

from collections import namedtuple
import sys

import numpy as np

from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.helpers import check_for_stop_sign


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


def dimer_method(geoms, calc_getter, N_init=None,
                 max_step=0.1, max_cycles=50,
                 trial_angle=5, angle_thresh=5, dR_base=0.01,
                 dx=0.001, alpha=0.1,
                 ana_2dpot=False):
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

        Default parameters from [1]
            max_step = 0.1 bohr
            dx = 0.01 bohr
            dR_base = = 0.01 bohr
    """
    # Parameters
    trial_rad = np.deg2rad(trial_angle)
    angle_thresh_rad = np.deg2rad(angle_thresh)

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
        N = N_init
        if N is None:
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
    forces0 = [geom0.forces, ]
    directions = []
    betas = []
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

        f0 = forces0[-1]
        # f0 = geom0.forces
        f0_rms = np.sqrt(np.power(f0, 2).mean())
        f0_max = f0.max()

        row_args = [i, C, f0_max, f0_rms]
        table.print_row(row_args)
        if f0_rms <= 1e-3 and f0_max <= 1.5e-3:
            table.print("Converged!")
            break

        # Translation
        f0_mod = get_f_mod(f0, N_rot, C)

        # Conjugate gradient
        if i == 0:
            directions.append(f0_mod/np.linalg.norm(f0_mod))

        direction = directions[-1]
        # Use constant alpha of 0.001
        alpha = dx
        coords0_trans = coords0 + alpha*direction
        geom0_star = geom_getter(coords0_trans)
        f0_mod_star = get_f_mod(geom0_star.forces, N_rot, C)
        forces0.append(f0_mod_star)
        beta = f0_mod_star.dot(f0_mod_star - f0_mod) / np.sqrt(f0_mod.dot(f0_mod))
        beta = max(beta, 0.0)
        print(f"beta {beta:.4f}")
        # betas.append(beta)
        direction = f0_mod_star + beta*direction
        direction /= np.linalg.norm(direction)
        directions.append(direction)

        # N0_mod = f0_mod / np.linalg.norm(f0_mod)
        # trial_coords0 = coords0 + dx*N0_mod
        # trial_geom0 = geom_getter(trial_coords0)
        # trial_f0 = trial_geom0.forces
        # trial_f0_mod = get_f_mod(trial_f0, N_rot, C)
        # F_mod = (trial_f0_mod + f0_mod).dot(N0_mod) / 2
        # C_mod = (trial_f0_mod - f0_mod).dot(N0_mod) / dx
        # step = (-F_mod/C_mod + dx/2)*N0_mod
        # step_norm = np.linalg.norm(step)
        # if step_norm > max_step:
            # step = max_step * f0_mod
            # table.print(f"Step norm {step_norm:.2f} is too large. Scaling down!")

        # # Small displacement of midpoint
        # coords0_trans = coords0 + step
        coords1_trans = coords0_trans + dR*N_rot
        coords2_trans = coords0_trans - dR*N_rot

        # Save cycle information
        org_coords = np.array((coords1, coords0, coords2))
        rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        trans_coords = np.array((coords1_trans, coords0_trans, coords2_trans))
        dc = DimerCycle(org_coords, trial_coords, rot_coords, trans_coords, f0, f0_mod)
        dimer_cycles.append(dc)

        # Update dimer coordinates for next cycle
        geom0.coords = coords0_trans
        geom1.coords = coords1_trans
        geom2.coords = coords2_trans
        N = N_rot

        header = "Cycle Curvature max(f0) rms(f0)".split()
        if ana_2dpot:
            hess = geom0.hessian
            w, v = np.linalg.eig(hess)
            trans_mode = v[:,0]
            mode_ovlp = trans_mode.dot(N_rot)
        if check_for_stop_sign():
            break
    return dimer_cycles


if __name__ == "__main__":
    run()
