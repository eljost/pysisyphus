from collections import namedtuple
import logging
import sys

import cloudpickle
import numpy as np

from pysisyphus.constants import EVANG2AUBOHR
from pysisyphus.helpers import check_for_stop_sign, get_geom_getter
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.optimizers.closures import lbfgs_closure_
import pysisyphus.optimizers.closures as closures
from pysisyphus.TablePrinter import TablePrinter


logger = logging.getLogger("tsoptimizer")


DimerCycle = namedtuple("DimerCycle",
                        "org_coords rot_coords trans_coords f0 f_tran",
)

DimerPickle = namedtuple("DimerPickle",
                         "coords0 N dR"
)

DimerResult = namedtuple("DimerResult",
                         "dimer_cycles force_evals geom0 converged")

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


def get_rms(arr):
    return np.sqrt(np.mean(np.square(arr)))


def get_lambda(f_parallel):
    rms = get_rms(f_parallel)

    if rms < 0.5*EVANG2AUBOHR:
        l_ = 1.
    elif rms < 1.0*EVANG2AUBOHR:
        l_ = 0.5
    elif rms < 2.0*EVANG2AUBOHR:
        l_ = 0.25
    else:
        l_ = 0.1
    return l_


def get_f_tran_mod(f, N, C):
    """Eq. (15) and (16) from [5]."""
    f_parallel = f.dot(N)*N
    f_perp = f - f_parallel
    if C < 0:
        lambda_ = get_lambda(f_parallel)
        f_tran = f_perp - lambda_ * f_parallel
    else:
        perp_rms =  get_rms(f_perp)
        if perp_rms < 2*EVANG2AUBOHR:
            f_tran = 0.5*f_perp - f_parallel
        else:
            f_tran = f_perp - 0.5*f_parallel
        # f_tran = -f_parallel
    return f_tran, f_parallel, f_perp


def get_f_tran_org(f, N, C):
    f_parallel = f.dot(N)*N
    f_perp = f - f_parallel
    if C < 0:
        f_tran = f_perp - f_parallel
    else:
        f_tran = -f_parallel
    return f_tran, f_parallel, f_perp

    # return f_tran
    # WTH; It seems to make a difference if I return f_tran
    # or f_mod here ... wtf? This only happens when MB is used
    # for the translation...
    # f_mod = -f.dot(N)*N
    # if C < 0:
        # f_mod = f + 2*f_mod
    # np.testing.assert_allclose(f_mod, f_tran)
    # print(f"diff={np.linalg.norm(f_tran-f_mod):.4e}")
    # print(f"type(f_mod)={type(f_mod)}")
    # return f_mod


def get_f_perp(f1, f2, N):
    f1_perp = perpendicular_force(f1, N)
    f2_perp = perpendicular_force(f2, N)
    f_perp = f1_perp - f2_perp
    return f_perp


def get_theta0(f1, f2, N):
    """Construct vector theta from f1 and f2."""
    f_perp = get_f_perp(f1, f2, N)
    theta = f_perp / np.linalg.norm(f_perp)
    return theta


def write_progress(geom0):
    ts_fn = "dimer_ts.xyz"
    ts_xyz_str = geom0.as_xyz()
    with open(ts_fn, "w") as handle:
        handle.write(ts_xyz_str)
    logger.debug(f"Wrote current TS geometry to '{ts_fn}'.")


def dimer_method(geoms, calc_getter, N_init=None,
                 max_step=0.1, max_cycles=50,
                 max_rots=10, interpolate=True,
                 rot_type="fourier",
                 rot_opt="lbfgs", trans_opt="lbfgs",
                 trans_memory=5,
                 trial_angle=5, angle_tol=0.5, dR_base=0.01,
                 restrict_step="scale", ana_2dpot=False,
                 f_thresh=1e-3, rot_f_thresh=2e-3,
                 zero_weights=[], dimer_pickle=None,
                 f_tran_mod=True,
                 multiple_translations=False, max_translations=10):
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
        # Modified Broyden
        [5] https://doi.org/10.1021/ct9005147

        To add:
            Comparing curvatures and adding π/2 if appropriate.

        Default parameters from [1]
            max_step = 0.1 bohr
            dR_base = = 0.01 bohr
    """
    # Parameters
    rad_tol = np.deg2rad(angle_tol)
    rms_f_thresh = float(f_thresh)
    max_f_thresh = 1.5 * rms_f_thresh
    max_translations = max_translations if multiple_translations else 1

    opt_name_dict = {
        "cg": "conjugate gradient",
        "lbfgs": "L-BFGS",
        "mb": "modified Broyden",
        "sd": "steepst descent",
    }
    assert rot_opt in opt_name_dict.keys()
    assert rot_type in "fourier direct".split()
    # Translation using L-BFGS as described in [4]
    # Modified broyden as proposed in [5] 10.1021/ct9005147
    trans_closures = {
        "lbfgs": closures.lbfgs_closure_,
        "mb": closures.modified_broyden_closure,
    }
    print(f"Using {opt_name_dict[rot_opt]} for dimer rotation.")
    print(f"Using {opt_name_dict[trans_opt]} for dimer translation.")
    print(f"Keeping information of last {trans_memory} cycles for "
           "translation optimization.")

    f_tran_dict = {
        True: get_f_tran_mod,
        False: get_f_tran_org
    }
    get_f_tran = f_tran_dict[f_tran_mod]

    rot_header = "Rot._cycle Curvature rms(f_rot)".split()
    rot_fmts = "int float float".split()
    rot_table = TablePrinter(rot_header, rot_fmts, width=16, shift_left=2)
    trans_header = "Trans._cycle Curvature max(|f0|) rms(f0) rms(step)".split()
    trans_fmts = "int float float float float".split()
    trans_table = TablePrinter(trans_header, trans_fmts, width=16)

    geom_getter = get_geom_getter(geoms[0], calc_getter)

    assert len(geoms) in (1, 2), "geoms argument must be of length 1 or 2!"
    # if dimer_pickle:
        # with open(dimer_pickle, "rb") as handle:
            # dimer_tuple = cloudpickle.load(handle)
    if len(geoms) == 2:
        geom1, geom2 = geoms
        dR = np.linalg.norm(geom1.coords - geom2.coords) / 2
        # N points from geom2 to geom1
        N = make_unit_vec(geom1.coords, geom2.coords)
        coords0 = geom2.coords + dR*N
        geom0 = geom_getter(coords0)
    # This handles cases where only one geometry is supplied. We use the
    # given geometry as midpoint, and use the given N_init. If N_init is
    # none, we select a random direction and derive geom1 and geom2 from it.
    else:
        geom0 = geoms[0]
        coords0 = geom0.coords
        # Assign random unit vector and use dR_base to for dR
        coord_size = geom0.coords.size
        if N_init is None:
            N_init  = np.random.rand(coord_size)
        N = np.array(N_init)
        if ana_2dpot:
            N[2] = 0
        N /= np.linalg.norm(N)
        coords1 = geom0.coords + dR_base*N
        geom1 = geom_getter(coords1)
        coords2 = geom0.coords - dR_base*N
        geom2 = geom_getter(coords2)
        dR = dR_base

    geom0.calculator.base_name += "_image0"
    geom1.calculator.base_name += "_image1"
    geom2.calculator.base_name += "_image2"

    dimer_pickle = DimerPickle(coords0, N, dR)
    with open("dimer_pickle", "wb") as handle:
        cloudpickle.dump(dimer_pickle, handle)

    dimer_cycles = list()

    print("Using N:", N)
    def f_tran_getter(coords, N, C):
        # The force for the given coord should be already set.
        np.testing.assert_allclose(geom0.coords, coords)
        # Only return f_tran and drop the parallel and perpendicular component.
        return get_f_tran(geom0.forces, N, C)[0]

    def rot_force_getter(N, f1, f2):
        return get_f_perp(f1, f2, N)

    def restrict_max_step_comp(x, step, max_step=max_step):
        step_max = np.abs(step).max()
        if step_max > max_step:
            factor = max_step / step_max
            step *= factor
        return step

    def restrict_step_length(x, step, max_step_length=max_step):
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

    tot_rot_force_evals = 0
    # Create the translation optimizer in the first cycle of the loop.
    trans_opt_kwargs = {
        "restrict_step": rstr_func,
        "M": trans_memory,
        # "beta": 0.01,
    }
    if trans_opt == "mb":
        trans_opt_kwargs["beta"] = 0.01
    trans_optimizer = trans_closures[trans_opt](f_tran_getter, **trans_opt_kwargs)

    def cbm_rot_force_getter(coords1, N, f1, f2):
        return get_f_perp(f1, f2, N)
    converged = False
    # In the first dimer cycle f0 and f1 aren't set and have to be calculated.
    # For cycles i > 0 f0 and f1 will be already set here.
    add_force_evals = 2
    for i in range(max_cycles):
        logger.debug(f"Dimer macro cycle {i:03d}")
        print(highlight_text(f"Dimer Cycle {i:03d}"))
        f0 = geom0.forces
        f1 = geom1.forces
        f2 = 2*f0 - f1

        coords0 = geom0.coords
        coords1 = geom1.coords
        coords2 = geom2.coords

        C = get_curvature(f1, f2, N, dR)

        f0_rms = get_rms(f0)
        f0_max = np.abs(f0).max()
        trans_table.print(f"@ {i:03d}: C={C:.6f}; max(|f0|)={f0_max:.6f}; rms(f0)={f0_rms:.6f}")
        print()
        converged = C < 0 and f0_rms <= rms_f_thresh and f0_max <= max_f_thresh
        if converged:
            rot_table.print("@ Converged!")
            break

        rot_force_evals = 0
        rot_force0 = get_f_perp(f1, f2, N)
        # Initialize some data structure for the rotation optimizers
        if rot_opt == "cg":
            # Lists for conjugate gradient
            G_perps = [rot_force0, ]
            F_rots = [rot_force0, ]
        # LBFGS optimizer
        elif rot_opt == "lbfgs" :
            rot_lbfgs = lbfgs_closure_(rot_force_getter)
        # Modified broyden
        elif rot_opt == "mb":
            rot_mb = closures.modified_broyden_closure(rot_force_getter)

        def cbm_restrict(coords1, step, coords0=coords0, dR=dR):
            """Constrain of R1 back on hypersphere."""
            coords1_rot = coords1 + step
            coords1_dir = coords1_rot - coords0
            coords1_dir /= np.linalg.norm(coords1_dir)
            coords1_rot = coords0 + coords1_dir*dR
            return coords1_rot - coords1
        rot_cbm_kwargs = {
            "restrict_step": cbm_restrict,
            "beta": 0.01,
            "M": 3,
        }
        rot_cbm = closures.modified_broyden_closure(cbm_rot_force_getter, **rot_cbm_kwargs)
        for j in range(max_rots):
            logger.debug(f"Rotation cycle {j:02d}")
            if rot_type == "fourier":
                C = get_curvature(f1, f2, N, dR)
                logger.debug(f"C={C:.6f}")
                # Theta from L-BFGS as implemented in DL-FIND dlf_dimer.f90
                if rot_opt == "lbfgs":
                    theta_dir, rot_force = rot_lbfgs(N, f1, f2)
                elif rot_opt == "mb":
                    theta_dir, rot_force = rot_mb(N, f1, f2)
                # Theta from conjugate gradient
                elif j > 0 and rot_opt == "cg":
                    # f_perp = weights.dot(get_f_perp(f1, f2, N))
                    f_perp = get_f_perp(f1, f2, N)
                    rot_force = f_perp
                    gamma = (f_perp - F_rots[-1]).dot(f_perp) / f_perp.dot(f_perp)
                    G_last = G_perps[-1]
                    # theta will be present with j > 0.
                    G_perp = f_perp + gamma * np.linalg.norm(G_last)*theta  # noqa: F821
                    theta_dir = G_perp
                    F_rots.append(f_perp)
                    G_perps.append(G_perp)
                # Theta from plain steepest descent F_rot/|F_rot|
                else:
                    # theta_dir = weights.dot(get_f_perp(f1, f2, N))
                    theta_dir = get_f_perp(f1, f2, N)
                # Remove component that is parallel to N
                theta_dir = theta_dir - theta_dir.dot(N)*N
                theta = theta_dir / np.linalg.norm(theta_dir)

                # Get rotated endpoint geometries. The rotation takes place in a plane
                # spanned by N and theta. Theta is a unit vector perpendicular to N that
                # can be formed from the perpendicular components of the forces at the
                # endpoints.

                # Derivative of the curvature, Eq. (29) in [2]
                # (f2 - f1) or -(f1 - f2)
                dC = 2*(f0-f1).dot(theta)/dR
                rad_trial = -0.5*np.arctan2(dC, 2*abs(C))
                logger.debug(f"rad_trial={rad_trial:.2f}")
                # print(f"rad_trial={rad_trial:.2f} rad_tol={rad_tol:.2f}")
                if np.abs(rad_trial) < rad_tol:
                    logger.debug(f"rad_trial={rad_trial:.2f} below threshold. Breaking.")
                    break

                # Trial rotation for finite difference calculation of rotational force
                # and rotational curvature.
                coords1_trial = rotate_R1(coords0, rad_trial, N, theta, dR)
                f1_trial = geom1.get_energy_and_forces_at(coords1_trial)["forces"]
                rot_force_evals += 1
                f2_trial = 2*f0 - f1_trial
                N_trial = make_unit_vec(coords1_trial, coords0)

                C_trial = get_curvature(f1_trial, f2_trial, N_trial, dR)

                b1 = 0.5 * dC
                a1 = (C - C_trial + b1*np.sin(2*rad_trial)) / (1-np.cos(2*rad_trial))
                a0 = 2 * (C - a1)

                rad_min = 0.5 * np.arctan(b1/a1)
                logger.debug(f"rad_min={rad_min:.2f}")
                def get_C(theta_rad):
                    return a0/2 + a1*np.cos(2*theta_rad) + b1*np.sin(2*theta_rad)
                C_min = get_C(rad_min)  # lgtm [py/multiple-definition]
                if C_min > C:
                    rad_min += np.deg2rad(90)
                    C_min_new = get_C(rad_min)
                    logger.debug( "Predicted theta_min lead us to a curvature maximum "
                                 f"(C(theta)={C_min:.6f}). Adding pi/2 to theta_min. "
                                 f"(C(theta+pi/2)={C_min_new:.6f})"
                    )
                    C_min = C_min_new

                # TODO: handle cases where the curvature is still positive, but
                # the angle is small, so the rotation is skipped.
                # Don't do rotation for small angles
                if np.abs(rad_min) < rad_tol:
                    logger.debug(f"rad_min={rad_min:.2f} below threshold. Breaking.")
                    break
                coords1_rot = rotate_R1(coords0, rad_min, N, theta, dR)

                # Interpolate force at coords1_rot; see Eq. (12) in [4]
                if interpolate:
                    f1 = (np.sin(rad_trial-rad_min)/np.sin(rad_trial)*f1
                           + np.sin(rad_min)/np.sin(rad_trial)*f1_trial
                           + (1 - np.cos(rad_min) - np.sin(rad_min)
                              * np.tan(rad_trial/2))*f0
                    )
                else:
                    f1 = geom1.get_energy_and_forces_at(coords1_rot)["forces"]
                    rot_force_evals += 1
            elif rot_type == "direct":
                rot_force = get_f_perp(f1, f2, N)
                rot_force_rms = get_rms(rot_force)
                # print(f"rot_force_rms={rot_force_rms:.4e}, thresh={rot_f_thresh:.4e}")
                if rot_force_rms <= rot_f_thresh:
                    break
                rot_step, rot_f = rot_cbm(coords1, N, f1, f2)
                coords1_rot = coords1 + rot_step
                geom1.coords = coords1_rot
                coords1 = coords1_rot
                f1 = geom1.forces
                rot_force_evals += 1
            else:
                raise NotImplementedError("Invalid 'rot_type'!")

            N = make_unit_vec(coords1_rot, coords0)
            coords2_rot = coords0 - N*dR
            f2 = 2*f0 - f1
            C = get_curvature(f1, f2, N, dR)
            rot_force = get_f_perp(f1, f2, N)
            rot_force_rms = get_rms(rot_force)
            logger.debug("")
            if j == 0:
                rot_table.print_header()
            rot_table.print_row((j, C, rot_force_rms))
        tot_rot_force_evals += rot_force_evals
        rot_str = (f"Did {rot_force_evals} force evaluation(s) and {j} "
                    "dimer rotation(s).")
        rot_table.print(rot_str)
        logger.debug(rot_str)
        print()

        # If multiple_translations == False then max_translations is 1
        # and we will do only one iteration.
        for trans_i in range(max_translations):
            _, f_parallel, f_perp = get_f_tran(f0, N, C)
            prev_f_par_rms = get_rms(f_parallel)
            # prev_f_perp_rms = get_rms(f_perp)
            prev_C = C
            f0_rms = get_rms(f0)
            f0_max = max(np.abs(f0))
            converged = C < 0 and f0_rms <= rms_f_thresh and f0_max <= max_f_thresh
            if converged:
                break
            step, f_tran = trans_optimizer(coords0, N, C)
            step_rms = get_rms(step)
            coords0_trans = coords0 + step
            coords1_trans = coords0_trans + dR*N
            coords2_trans = coords0_trans - dR*N
            geom0.coords = coords0_trans
            geom1.coords = coords1_trans
            geom2.coords = coords2_trans
            coords0 = geom0.coords
            # Calculate new forces for translated dimer
            f0 = geom0.forces
            f1 = geom1.forces
            add_force_evals += 2
            f2 = 2*f0 - f1
            C = get_curvature(f1, f2, N, dR)
            f0_rms = get_rms(f0)
            f0_max = max(np.abs(f0))
            if multiple_translations:
                if trans_i == 0:
                    trans_table.print_header()
                trans_args = (trans_i, C, f0_max, f0_rms, step_rms)
                trans_table.print_row(trans_args)
            else:
                trans_table.print("Did dimer translation.")
            _, f_parallel, f_perp = get_f_tran(f0, N, C)
            f_par_rms = get_rms(f_parallel)
            # f_perp_rms = get_rms(f_perp)
            # Check for sign change of curvature
            if (prev_C < 0) and (np.sign(C/prev_C) < 0):
                trans_table.print("Curvature became positive!")
                break
            if (C < 0) and (f_par_rms > prev_f_par_rms):
                break
            # elif ((C > 0) and (f_par_rms < prev_f_par_rms)
                  # or (f_perp_rms > prev_f_perp_rms)):
            elif (C > 0) and (f_par_rms < prev_f_par_rms):
                break
            prev_f_par_rms = f_par_rms  # lgtm [py/multiple-definition]
            # prev_f_perp_rms = f_perp_rms

        # Save cycle information
        org_coords = np.array((coords1, coords0, coords2))
        try:
            rot_coords = np.array((coords1_rot, coords0, coords2_rot))
        except NameError:
            rot_coords = np.array((coords1, coords0, coords2))
        trans_coords = np.array((coords1_trans, coords0_trans, coords2_trans))
        dc = DimerCycle(org_coords, rot_coords, trans_coords, f0, f_tran)
        dimer_cycles.append(dc)

        write_progress(geom0)

        if check_for_stop_sign():
            break
        logger.debug("")
        print()
        sys.stdout.flush()
    print(f"@Did {tot_rot_force_evals} force evaluations in the rotation steps "
          f"using the {opt_name_dict[rot_opt]} optimizer.")
    tot_force_evals = tot_rot_force_evals + add_force_evals
    print(f"@Used {add_force_evals} additional force evaluations for a total of "
          f"{tot_force_evals} force evaluations.")

    dimer_results = DimerResult(dimer_cycles=dimer_cycles,
                                force_evals=tot_force_evals,
                                geom0=geom0,
                                converged=converged)

    return dimer_results
