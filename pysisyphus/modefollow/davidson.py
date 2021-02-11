# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1523908
#     Neugebauer, Reiher 2002
# [2] https://reiher.ethz.ch/software/akira.html


from collections import namedtuple
import sys

import numpy as np

from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.Geometry import get_trans_rot_projector
from pysisyphus.modefollow.NormalMode import NormalMode
from pysisyphus.TablePrinter import TablePrinter


DavidsonResult = namedtuple(
    "DavidsonResult",
    "cur_cycle converged final_modes qs nus mode_inds res_rms",
)


def forces_fin_diff(forces_getter, coords, b, step_size):
    plus = forces_getter(coords + b)
    minus = forces_getter(coords - b)
    fd = (minus - plus) / (2 * step_size)
    return fd


def block_davidson(
    cart_coords,
    masses,
    forces_getter,
    guess_modes,
    lowest=None,
    trial_step_size=0.01,
    hessian_precon=None,
    max_cycles=25,
    res_rms_thresh=1e-4,
    start_precon=5,
    remove_trans_rot=True,
    print_level=1,
):
    num = len(guess_modes)
    B_full = np.zeros((len(guess_modes[0]), num * max_cycles))
    S_full = np.zeros_like(B_full)
    I = np.eye(cart_coords.size)
    masses_rep = np.repeat(masses, 3)
    msqrt = np.sqrt(masses_rep)

    # Projector to remove translation and rotation
    P = get_trans_rot_projector(cart_coords, masses)
    guess_modes = [
        NormalMode(P.dot(mode.l_mw) / msqrt, masses_rep) for mode in guess_modes
    ]

    col_fmts = "int int int float_short float str".split()
    header = ("#", "subspace size", "mode", "ṽ / cm⁻¹", "rms(r)", "Conv")
    fmts_update = {"float_short": "{: >11.2f}"}
    table = TablePrinter(header, col_fmts, width=11, fmts_update=fmts_update)
    if print_level == 1:
        table.print_header()

    b_prev = np.array([mode.l_mw for mode in guess_modes]).T
    for i in range(max_cycles):
        # Add new basis vectors to B matrix
        b = np.array([mode.l_mw for mode in guess_modes]).T
        from_ = i * num
        to_ = (i + 1) * num
        B_full[:, from_:to_] = b

        # Estimate action of Hessian by finite differences.
        for j in range(num):
            mode = guess_modes[j]
            # Get a step size in mass-weighted coordinates that results
            # in the desired 'trial_step_size' in not-mass-weighted coordinates.
            mw_step_size = mode.mw_norm_for_norm(trial_step_size)
            # Actual step in non-mass-weighted coordinates
            step = trial_step_size * mode.l
            S_full[:, from_ + j] = (
                # Convert to mass-weighted coordinates
                forces_fin_diff(forces_getter, cart_coords, step, mw_step_size)
                / msqrt
            )

        # Views on columns that are actually set
        B = B_full[:, :to_]
        S = S_full[:, :to_]

        # Calculate and symmetrize approximate hessian
        Hm = B.T.dot(S)
        Hm = (Hm + Hm.T) / 2
        # Diagonalize small Hessian
        w, v = np.linalg.eigh(Hm)

        # Approximations to exact eigenvectors in current cycle
        approx_modes = (v * B[:, :, None]).sum(axis=1).T
        # Calculate overlaps between previous root and the new approximate
        # normal modes for root following.
        if lowest is None:
            # 2D overlap array. approx_modes in row, b_prev in columns.
            overlaps = np.einsum("ij,jk->ik", approx_modes, b_prev)
            mode_inds = np.abs(overlaps).argmax(axis=0)
        else:
            mode_inds = np.arange(lowest)
        b_prev = approx_modes[mode_inds].T

        # Eq. (7) in [1]
        residues = (v * (S[:, :, None] - w * B[:, :, None])).sum(axis=1)

        # Determine preconditioner matrix
        #
        # Use supplied matrix
        if hessian_precon is not None:
            precon_mat = hessian_precon
        # Reconstruct Hessian, but only start after some cycles
        elif i >= start_precon:
            precon_mat = B.dot(Hm).dot(B.T)
        # No preconditioning if no matrix was supplied or we are in an early cycle.
        else:
            precon_mat = None

        # Construct new basis vector from residuum of selected mode
        b = np.zeros_like(b_prev)
        for j, mode_ind in enumerate(mode_inds):
            r = residues[:, mode_ind]
            if precon_mat is not None:
                # Construct actual preconditioner X
                X = np.linalg.pinv(precon_mat - w[mode_ind] * I, rcond=1e-8)
                b[:, j] = X.dot(r)
            else:
                b[:, j] = r

        # Project out translation and rotation from new mode guess
        if remove_trans_rot:
            b = P.dot(b)
        # Orthogonalize new vectors against preset vectors
        b = np.linalg.qr(np.concatenate((B, b), axis=1))[0][:, -num:]

        # New NormalMode from non-mass-weighted displacements
        guess_modes = [NormalMode(b_ / msqrt, masses_rep) for b_ in b.T]

        # Calculate wavenumbers
        nus = eigval_to_wavenumber(w)

        # Check convergence criteria
        max_res = np.abs(residues).max(axis=0)
        res_rms = np.sqrt(np.mean(residues ** 2, axis=0))

        converged = res_rms < res_rms_thresh
        # Print progress if requested
        if print_level == 2:
            print(f"Cycle {i:02d}")
            print("\t #  |    ṽ / cm⁻¹|   rms(r)   | max(|r|) ")
            print("\t------------------------------------------")
            for j, (nu, rms, mr) in enumerate(zip(nus, res_rms, max_res)):
                sel_str = "*" if (j in mode_inds) else " "
                conv_str = "✓" if converged[j] else ""
                print(
                    f"\t{j:02d}{sel_str} | {nu:> 10.2f} | {rms:.8f} | {mr:.8f} {conv_str}"
                )
            print()
        elif print_level == 1:
            for j in mode_inds:
                conv_str = "✓" if converged[j] else "✗"
                table.print_row((i, B.shape[1], j, nus[j], res_rms[j], conv_str))

        # Convergence is signalled using only the roots we are actually interested in
        modes_converged = all(converged[mode_inds])
        if modes_converged:
            if print_level > 0:
                print(f"\tDavidson procedure converged in {i+1} cycles!")
                if lowest is not None:
                    nus_str = np.array2string(nus[mode_inds], precision=2)
                    print(f"\tLowest {lowest} wavenumbers: {nus_str} cm⁻¹")
                    neg_nus = sum(nus[mode_inds] < 0)
                    type_ = "minimum" if (neg_nus == 0) else f"saddle point of index {neg_nus}"
                    print(f"\tThis geometry seems to be a {type_} on the PES.")
            break
        sys.stdout.flush()

    result = DavidsonResult(
        cur_cycle=i,
        converged=modes_converged,
        final_modes=guess_modes,
        qs=approx_modes,
        nus=nus,
        mode_inds=mode_inds,
        res_rms=res_rms,
    )

    return result


def geom_davidson(geom, *args, **kwargs):
    def forces_getter(cart_coords):
        return geom.get_energy_and_cart_forces_at(cart_coords)["forces"]

    return block_davidson(geom.cart_coords, geom.masses, forces_getter, *args, **kwargs)
