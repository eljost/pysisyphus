#!/usr/bin/env python3

from collections import namedtuple

import numpy as np

from pysisyphus.helpers import check_for_stop_sign, get_geom_getter
from pysisyphus.optimizers.closures import lbfgs_closure
from pysisyphus.TablePrinter import TablePrinter


DimerCycle = namedtuple("DimerCycle",
                        "org_coords rot_coords trans_coords f0 f0_mod",
)

def make_unit_vec(vec1, vec2):
    """Return unit vector pointing from vec2 to vec1."""
    diff = vec1 - vec2
    return diff / np.linalg.norm(diff)


def perpendicular_force(force, vec):
    """Return the perpendicular component of force along vec."""
    return force - force.dot(vec)*vec


def rotate_R1(coords0, rad, N, theta, dR):
    """Rotate dimer at coords0 to get a new coords1.

    Only valid for rotation of R1!"""
    return coords0 + (N*np.cos(rad) + theta*np.sin(rad)) * dR


def get_curvature(f1, f2, N, dR):
    return (f2 - f1).dot(N) / (2*dR)


def get_f_mod(f, N, C):
    f_mod = -f.dot(N)*N
    if C < 0:
        f_mod = f + 2*f_mod
    return f_mod


def make_theta(f1, f2, N):
    """Construct vector theta from f1 and f2."""
    f1_perp = perpendicular_force(f1, N)
    f2_perp = perpendicular_force(f2, N)
    f_perp = f1_perp - f2_perp
    theta = f_perp / np.linalg.norm(f_perp)
    return theta


def write_progress(geom0):
    ts_fn = "dimer_ts.xyz"
    ts_xyz_str = geom0.as_xyz()
    with open(ts_fn, "w") as handle:
        handle.write(ts_xyz_str)
    print(f"Wrote current TS geometry to '{ts_fn}'.")


def dimer_method(geoms, calc_getter, N_init=None,
                 max_step=0.1, max_cycles=50,
                 trial_angle=5, angle_thresh=0.5, dR_base=0.01,
                 restrict_step="scale", ana_2dpot=False,
                 f_thresh=1e-3):
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
            dR_base = = 0.01 bohr
    """
    # Parameters
    angle_thresh_rad = np.deg2rad(angle_thresh)
    rms_f_thresh = float(f_thresh)
    max_f_thresh = 1.5 * rms_f_thresh

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
        coords0 = geom0.coords
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

    dimer_cycles = list()
    coords0_list = [coords0.copy(), ]
    N_list = [N.copy(), ]

    print("Using N:", N)
    def f0_mod_getter(coords, N, C):
        # results = geom0.get_energy_and_forces_at(coords)
        # forces = results["forces"]
        # return get_f_mod(forces, N, C)
        geom0.coords = coords
        forces = geom0.forces
        return get_f_mod(forces, N, C)

    def restrict_max_step_comp(step, max_step=max_step):
        step_max = np.abs(step).max()
        if step_max > max_step:
            factor = max_step / step_max
            step *= factor
        return step

    def restrict_step_length(step, max_step_length=max_step):
        step_norm = np.linalg.norm(step)
        if step_norm > max_step_length:
            step_direction = step / step_norm
            step = step_direction * max_step_length
        return step

    rstr_dict = {
        "max": restrict_max_step_comp,
        "scale": restrict_step_length,
    }
    rstr_func = rstr_dict[restrict_step]

    table.print_header()
    for i in range(max_cycles):
        f0 = geom0.forces
        f1 = geom1.forces
        f2 = 2*f0 - f1

        coords0 = geom0.coords
        coords1 = geom1.coords
        coords2 = geom2.coords

        C = get_curvature(f1, f2, N, dR)

        f0_rms = np.sqrt(np.power(f0, 2).mean())
        f0_max = np.abs(f0).max()
        row_args = [i, C, f0_max, f0_rms]
        table.print_row(row_args)
        if C < 0 and f0_rms <= rms_f_thresh and f0_max <= max_f_thresh:
            write_progress(geom0)
            table.print("Converged!")
            break

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.
        theta = make_theta(f1, f2, N)

        # Derivative of the curvature, Eq. (29) in [2]
        # (f2 - f1) or -(f1 - f2)
        dC = 2*(f0-f1).dot(theta)/dR
        rad_trial = -0.5*np.arctan2(dC, 2*abs(C))

        # Trial rotation for finite difference calculation of rotational force
        # and rotational curvature.
        coords1_trial = rotate_R1(coords0, rad_trial, N, theta, dR)
        f1_trial = geom1.get_energy_and_forces_at(coords1_trial)["forces"]
        f2_trial = 2*f0 - f1_trial
        N_trial = make_unit_vec(coords1_trial, coords0)
        coords2_trial = coords0 - N_trial*dR

        C_trial = get_curvature(f1_trial, f2_trial, N_trial, dR)
        theta_trial = make_theta(f1_trial, f2_trial, N_trial)

        b1 = 0.5 * dC
        a1 = (C - C_trial + b1*np.sin(2*rad_trial)) / (1-np.cos(2*rad_trial))
        a0 = 2 * (C - a1)

        rad_min = 0.5 * np.arctan(b1/a1)
        C_min = a0/2 + a1*np.cos(2*rad_min) + b1*np.sin(2*rad_min)
        if C_min > C:
            rad_min += np.deg2rad(90)

        # TODO: Multiple rotations after another with L-BFGS
        f1_extrapol = (np.sin(rad_trial-rad_min)/np.sin(rad_trial)*f1
                       + np.sin(rad_min)/np.sin(rad_trial)*f1_trial
                       + (1 - np.cos(rad_min) - np.sin(rad_min)
                          * np.tan(rad_trial/2))*f0
        )

        # TODO: handle cases where the curvature is still positive, but
        # the angle is small, so the rotation is skipped.
        if np.abs(rad_min) > angle_thresh_rad:
            coords1_rot = rotate_R1(coords0, rad_min, N, theta, dR)
            N = make_unit_vec(coords1_rot, coords0)
            coords2_rot = coords0 - N*dR
        # Don't do rotation for small angles
        else:
            table.print("Rotation angle too small. Skipping rotation.")

        f0_mod = get_f_mod(f0, N, C_min)

        # Create the translation optimizer in the first cycle of the loop.
        if i == 0:
            trans_lbfgs = lbfgs_closure(f0_mod, f0_mod_getter,
                                        restrict_step=rstr_func)

        if C > 0:
            N_trans = f0_mod.copy()
            N_trans /= np.linalg.norm(N_trans)
            step = max_step*N_trans
            # TODO: geom0 coords aren't updated here ...
            coords0_trans = coords0 + step
            geom0.coords = coords0_trans
        else:
            # Translation using L-BFGS as described in [4]
            coords0_trans, step, f0 = trans_lbfgs(coords0, N, C_min)

        # The coordinates of geom0 were already updated in the f0_mod_getter
        # call.
        coords1_trans = coords0_trans + dR*N
        coords2_trans = coords0_trans - dR*N

        # Save cycle information
        org_coords = np.array((coords1, coords0, coords2))
        rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        trans_coords = np.array((coords1_trans, coords0_trans, coords2_trans))
        dc = DimerCycle(org_coords, rot_coords, trans_coords, f0, f0_mod)
        dimer_cycles.append(dc)

        # Update dimer coordinates for next cycle
        # geom0.coords = coords0_trans
        geom1.coords = coords1_trans
        geom2.coords = coords2_trans

        if check_for_stop_sign():
            write_progress(geom0)
            break

    return dimer_cycles


if __name__ == "__main__":
    run()
