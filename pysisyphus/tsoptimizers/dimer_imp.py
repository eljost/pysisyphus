#!/usr/bin/env python3

from collections import namedtuple
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.helpers import check_for_stop_sign

np.set_printoptions(suppress=True, precision=4)


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


def cg_closure(first_force):
    prev_force = first_force
    prev_p = prev_force / np.linalg.norm(prev_force)

    def cg(cur_force):
        nonlocal prev_force
        nonlocal prev_p
        # Polak-Ribiere formula
        beta = cur_force.dot(cur_force - prev_force) / cur_force.dot(cur_force)
        beta = max(beta, 0)
        cur_p = cur_force + beta*prev_p

        prev_force = prev_force
        prev_p = cur_p

        return cur_p
    return cg


def bfgs_multiply(s_list, y_list, force):
    """Get a L-BFGS step.
    
    Algorithm 7.4 Nocedal, Num. Opt., p. 178."""
    q = -force
    cycles = len(s_list)
    alphas = list()
    rhos = list()
    # Store rho and alphas as they are also needed in the second loop
    for i in reversed(range(cycles)):
        s = s_list[i]
        y = y_list[i]
        rho = 1/y.dot(s)
        rhos.append(rho)
        alpha = rho * s.dot(q)
        alphas.append(alpha)
        q = q - alpha*y
    # Restore original order, so that rho[i] = 1/s_list[i].dot(y_list[i]) etc.
    alphas = alphas[::-1]
    rhos = rhos[::-1]

    r = q
    if cycles > 0:
        s = s_list[-1]
        y = y_list[-1]
        gamma = s.dot(y) / y.dot(y)
        r = gamma * q
    for i in range(cycles):
        s = s_list[i]
        y = y_list[i]
        beta = rhos[i] * y.dot(r)
        r = r + s*(alphas[i] - beta)

    return r


def lbfgs_closure(first_force, force_getter, m=10, restrict_step=None):
    s_list = list()
    y_list = list()
    forces = [first_force, ]
    cur_cycle = 0

    if restrict_step is None:
        restrict_step = lambda x: x

    def lbfgs(x, *getter_args):
        nonlocal cur_cycle
        nonlocal s_list
        nonlocal y_list

        prev_forces = forces[-1]
        step = -bfgs_multiply(s_list, y_list, prev_forces)
        step = restrict_step(step)
        new_x = x + step
        new_forces = force_getter(new_x, *getter_args)
        s = new_x - x
        s_list.append(s)
        y = prev_forces - new_forces
        y_list.append(y)
        forces.append(new_forces)
        if cur_cycle > m:
            s_list = s_list[:-m]
            y_list = y_list[:-m]
        cur_cycle += 1
        return new_x, step, new_forces
    return lbfgs


def make_theta(f1, f2, N):
    f1_perp = perpendicular_force(f1, N)
    f2_perp = perpendicular_force(f2, N)
    f_perp = f1_perp - f2_perp
    theta = f_perp / np.linalg.norm(f_perp)
    return theta


def plot_modes(trans_mode, N, N_trial, N_rot, coords0):
    fig, ax = plt.subplots()

    x, y = coords0[:2]
    trans_x, trans_y = trans_mode[:2]
    xt = np.array((x, x+trans_x))
    yt = np.array((y, y+trans_y))
    trans_line = mpl.lines.Line2D(xdata=xt, ydata=yt, label="Trans mode")
    lbls = "N N_trial N_rot".split()
    for i, (N_, l) in enumerate(zip((N, N_trial, N_rot), lbls)):
        xn, yn = N_[:2]
        xn_data = np.array((x, x+xn))
        yn_data = np.array((y, y+yn))
        n_line = mpl.lines.Line2D(xn_data, yn_data, label=l, color=f"C{i}")
        ax.add_line(n_line)

    # dimer_x, dimer_y = N[:2]
    # xn = np.array((x, x+dimer_x))
    # yn = np.array((y, y+dimer_y))
    # dimer_line = mpl.lines.Line2D(xdata=xn, ydata=yn, label="Dimer", color="red")

    # trial_x, trial_y = N_trial[:2]
    # xt_ = np.array((x, x+trial_x))
    # yt_ = np.array((y, y+trial_y))
    # trial_line = mpl.lines.Line2D(xdata=xt_, ydata=yt_, label="Trial", color="green")

    ax.scatter(x, y)
    ax.add_line(trans_line)
    # ax.add_line(dimer_line)
    # ax.add_line(trial_line)
    ax.legend()
    ax.set_xlim(-0.5, 2)
    ax.set_ylim(0, 3)

    plt.show()


def dimer_method(geoms, calc_getter, N_init=None,
                 max_step=0.1, max_cycles=50,
                 trial_angle=5, angle_thresh=0.5, dR_base=0.01,
                 dx=0.01, alpha=0.1,
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
    # trial_rad = np.deg2rad(trial_angle)
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
    directions = []
    betas = []
    def f0_mod_getter(coords, N, C):
        results = geom0.get_energy_and_forces_at(coords)
        forces = results["forces"]
        return get_f_mod(forces, N, C)

    def restrict_step(step, max_step=max_step):
        step_max = np.abs(step).max()
        if step_max > max_step:
            factor = max_step / step_max
            step *= factor
        return step

    for i in range(max_cycles):
        f0 = geom0.forces
        f1 = geom1.forces
        f2 = 2*f0 - f1

        coords0 = geom0.coords
        coords1 = geom1.coords
        coords2 = geom2.coords

        # Get rotated endpoint geometries. The rotation takes place in a plane
        # spanned by N and theta. Theta is a unit vector perpendicular to N that
        # can be formed from the perpendicular components of the forces at the
        # endpoints.
        theta = make_theta(f1, f2, N)

        C = get_curvature(f1, f2, N, dR)
        C_ = (f0 - f1).dot(N)/dR
        assert C == C_
        # Derivative of the curvature, Eq. (29) in [2]
        # (f2 - f1) or -(f1 - f2)
        dC = (f2 - f1).dot(theta) / dR
        dC_ = 2*(f0-f1).dot(theta)/dR
        assert dC == dC_
        # table.print(f"C={C:.4f}, dC={dC:.4f}")
        trial_rad = -0.5*np.arctan2(dC_, 2*abs(C))
        # print("trial_rad", np.rad2deg(trial_rad), "°")

        if ana_2dpot:
            hess = geom0.hessian
            w, v = np.linalg.eig(hess)
            trans_mode = v[:,0]
            mode_ovlp = trans_mode.dot(N)
            # print("mode_ovlp", mode_ovlp)

        # Trial rotation for finite difference calculation of rotational force
        # and rotational curvature.
        coords1_star = rotate_R1(coords0, trial_rad, N, theta, dR)
        f1_star = geom1.get_energy_and_forces_at(coords1_star)["forces"]
        f2_star = 2*f0 - f1_star
        N_star = make_unit_vec(coords1_star, coords0)
        coords2_star = coords0 - N_star*dR

        trial_coords = np.array((coords1_star, coords0, coords2_star))
        C_star = get_curvature(f1_star, f2_star, N_star, dR)
        # print("C_star", C_star)
        theta_star = make_theta(f1_star, f2_star, N_star)
        # dC_star = (f2_star - f1_star).dot(theta_star) / dR

        # a1 = (dC * np.cos(2*trial_rad) - dC_star) / (2*np.sin(2*trial_rad))
        b1 = 0.5 * dC
        a1 = (C - C_star + b1*np.sin(2*trial_rad)) / (1-np.cos(2*trial_rad))
        a0 = 2 * (C - a1)

        rad_min = 0.5 * np.arctan(b1/a1)
        # print("rad_min", np.rad2deg(rad_min), "°")
        C_min = a0/2 + a1*np.cos(2*rad_min) + b1*np.sin(2*rad_min)
        # print("C_min", C_min)
        if C_min > C:
            rad_min += np.deg2rad(90)
        # table.print("C_min", C_min)
        # table.print(f"trial: {np.rad2deg(trial_rad):.1f}°, rot: {np.rad2deg(rad_min):.1f}°")

        if np.abs(rad_min) > angle_thresh_rad:
            coords1_rot = rotate_R1(coords0, rad_min, N, theta, dR)
            N_rot = make_unit_vec(coords1_rot, coords0)
            coords2_rot = coords0 - N_rot*dR
            geom1.coords = coords1_rot
            geom2.coords = coords2_rot
        # Don't do rotation for small angles
        else:
            N_rot = N
            coords1_rot = coords1
            coords2_rot = coords2
            table.print("Rotation angle too small. Skipping rotation.")

        # plot_modes(trans_mode, N, N_star, N_rot, coords0)

        # Translation
        f0_mod = get_f_mod(f0, N_rot, C_min)
        # print("f0_mod", f0_mod)

        # Initialize conjugate gradient optimizer
        if i == 0:
            cg = cg_closure(f0_mod)
            trans_lbfgs = lbfgs_closure(f0_mod, f0_mod_getter,
                                        restrict_step=restrict_step)
        # Use conjugate gradient to determine the step direction
        # direction = cg(f0_mod)
        N_trans = f0_mod.copy()
        N_trans /= np.linalg.norm(N_trans)
        # import pdb; pdb.set_trace()

        if C > 0:
            step = max_step*N_trans
            table.print("curv positive!")
        else:
            # trial_coords0 = coords0 + dx*N_trans
            # trial_f0 = geom0.get_energy_and_forces_at(trial_coords0)["forces"]
            # trial_f0_mod = get_f_mod(trial_f0, N_rot, C_min)
            # F_mod = (trial_f0_mod + f0_mod).dot(N_trans) / 2
            # C_mod = (trial_f0_mod - f0_mod).dot(N_trans) / dx
            # step_lengths = (-F_mod/C_mod)# + dx/2)*N0_trans
            # table.print(f"step_lengths {step_lengths}")
            # step = step_lengths * N_trans
            # print("step", step)
            # if np.linalg.norm(step) > max_step:
                # step = N_trans * max_step
                # table.print("scaling down step")
            bfgs_x, bfgs_step, bfgs_force = trans_lbfgs(coords0, N_rot, C_min)
            # print()
            # print("BFGS")
            # print(bfgs_step)
            # print()
            step = bfgs_step
        # max_step_comp = np.abs(step).max()
        # table.print(f"max step comp is {max_step_comp}")
        # if max_step_comp > max_step:
            # factor = max_step / max_step_comp
            # step *= factor
            # table.print(f"scaled down with factor {factor}")

        # # Small displacement of midpoint
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
        geom0.coords = coords0_trans
        geom1.coords = coords1_trans
        geom2.coords = coords2_trans
        N = N_rot

        f0_rms = np.sqrt(np.power(f0, 2).mean())
        f0_max = np.abs(f0).max()
        row_args = [i, C, f0_max, f0_rms]
        table.print_row(row_args)
        if f0_rms <= 1e-3 and f0_max <= 1.5e-3:
            table.print("Converged!")
            break

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
