import numpy as np

from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.eval import eval_primitives
from pysisyphus.intcoords.exceptions import NeedNewInternalsException
from pysisyphus.intcoords.valid import dihedral_valid


def correct_dihedrals(new_dihedrals, old_dihedrals):
    """Dihedrals are periodic. Going from -179° to 179° is not a step of 358°,
    but a step of 2°. By considering the actual distance of the dihedrals from
    π the correct step can be calculated.

        dihedral step length  = abs(abs(new_dihedral) - π) + abs(abs(old_dihedral) - π)

    or put differently

        dihedral step length = abs(abs(new_dihedral - old_dihedral) - 2*π)

    The sign is left to be determined. Going from -179° to 179° (roughly π - -π = 2π)
    is a counter clockwise rotation and the dihedral has to decrease below -π. Going
    from 179° to -179° (roughly -π - π = -2π) is a clockwise rotation and the dihedral
    increases abvove π. So the correct sign corresponds to the negative sign of the
    original difference.

    original difference  2π -> dihedral must decrease -> sign = -1
    original difference -2π -> dihedral must increase -> sign = +1

    Overall, the old dihedral is modified by the actual step length with the correct
    sign."""
    new_dihedrals = np.atleast_1d(new_dihedrals)
    old_dihedrals = np.atleast_1d(old_dihedrals)
    dihedrals_step = new_dihedrals - old_dihedrals
    shifted_by_2pi = np.abs(np.abs(dihedrals_step) - 2 * np.pi) < np.pi / 2
    corrected_dihedrals = new_dihedrals.copy()
    corrected_dihedrals[shifted_by_2pi] -= (
        2 * np.pi * np.sign(dihedrals_step[shifted_by_2pi])
    )
    return corrected_dihedrals


def update_internals(
    new_coords3d,
    old_internals,
    primitives,
    dihedral_inds,
    rotation_inds,
    check_dihedrals=False,
    rotation_thresh=0.9,
    logger=None,
):
    prim_internals = eval_primitives(new_coords3d, primitives)
    new_internals = np.array([prim_int.val for prim_int in prim_internals])
    internal_diffs = new_internals - old_internals

    new_rotations = new_internals[rotation_inds]
    # Check for approaching singularity as discussed in the geomeTRIC paper. The
    # original code seems to check this only for linear molecules and instead
    # does some +2π/-2π magic, similar to how dihedrals differences are handled.
    if (np.abs(new_rotations / np.pi) >= rotation_thresh).any():
        raise NeedNewInternalsException(new_coords3d)

    dihedrals = [prim_internals[i] for i in dihedral_inds]
    if len(dihedrals) == 0:
        return prim_internals

    dihedral_diffs = internal_diffs[dihedral_inds]

    # Find differences that are shifted by 2*pi
    shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2 * np.pi) < np.pi / 2
    new_dihedrals = np.array([dihed.val for dihed in dihedrals])
    if any(shifted_by_2pi):
        new_dihedrals[shifted_by_2pi] -= (
            2 * np.pi * np.sign(dihedral_diffs[shifted_by_2pi])
        )
        # Update values
        for dihed, new_val in zip(dihedrals, new_dihedrals):
            dihed.val = new_val
    # See if dihedrals became invalid (collinear atoms)
    if check_dihedrals:
        are_valid = [dihedral_valid(new_coords3d, prim.inds) for prim in dihedrals]
        try:
            first_dihedral = dihedral_inds[0]
        except IndexError:
            first_dihedral = 0
        invalid_inds = [
            i + first_dihedral for i, is_valid in enumerate(are_valid) if not is_valid
        ]
        if len(invalid_inds) > 0:
            invalid_prims = [primitives[i] for i in invalid_inds]
            log(logger, "Dihedral(s) became invalid! Need new internal coordinates!")
            raise NeedNewInternalsException(
                new_coords3d, invalid_inds=invalid_inds, invalid_prims=invalid_prims
            )

    return prim_internals


def transform_int_step(
    int_step,
    old_cart_coords,
    cur_internals,
    Bt_inv_prim,
    primitives,
    dihedral_inds,
    rotation_inds,
    check_dihedrals=False,
    freeze_atoms=None,
    cart_rms_thresh=1e-6,
    Bt_inv_prim_getter=None,
    logger=None,
):
    """Transformation is done in primitive internals, so int_step must be given
    in primitive internals and not in DLC!"""

    if freeze_atoms is None:
        freeze_atoms = list()
    freeze_atoms = np.array(freeze_atoms, dtype=int)
    new_cart_coords = old_cart_coords.copy()
    remaining_int_step = int_step
    target_internals = cur_internals + int_step

    last_rms = 9999
    old_internals = cur_internals
    backtransform_failed = True
    for i in range(25):
        if Bt_inv_prim_getter is not None:
            Bt_inv_prim = Bt_inv_prim_getter(new_cart_coords)
            log(logger, f"Recalculated (B^T)^+ in microcycle {i}")

        cart_step = Bt_inv_prim.T.dot(remaining_int_step)
        # Remove step from frozen atoms.
        cart_step.reshape(-1, 3)[freeze_atoms] = 0.0
        cart_rms = np.sqrt(np.mean(cart_step ** 2))
        # Update cartesian coordinates
        new_cart_coords += cart_step
        # Determine new internal coordinates
        new_prim_ints = update_internals(
            new_cart_coords.reshape(-1, 3),
            old_internals,
            primitives,
            dihedral_inds,
            rotation_inds,
            check_dihedrals=check_dihedrals,
            logger=logger,
        )
        new_internals = [prim.val for prim in new_prim_ints]
        remaining_int_step = target_internals - new_internals
        internal_rms = np.sqrt(np.mean(remaining_int_step ** 2))
        log(
            logger,
            f"Cycle {i}: rms(Δcart)={cart_rms:1.4e}, rms(Δint.) = {internal_rms:1.5e}",
        )

        # This assumes the first cart_rms won't be > 9999 ;)
        if cart_rms < last_rms:
            # Store results of the conversion cycle for laster use, if
            # the internal-cartesian-transformation goes bad.
            best_cycle = (new_cart_coords.copy(), new_internals.copy())
            best_cycle_ind = i
        elif i != 0:
            # If the conversion somehow fails we fallback to the best previous step.
            log(logger, f"Backconversion failed! Falling back to step {best_cycle_ind}")
            new_cart_coords, new_internals = best_cycle
            break
        else:
            raise Exception(
                "Internal-cartesian back-transformation already "
                "failed in the first step. Aborting!"
            )
        old_internals = new_internals

        last_rms = cart_rms
        if cart_rms < cart_rms_thresh:
            log(
                logger,
                f"Internal->Cartesian transformation converged in {i+1} cycle(s)!",
            )
            backtransform_failed = False
            break

    # if check_dihedrals and (
    # not dihedrals_are_valid(new_cart_coords.reshape(-1, 3), dihedral_inds)
    # ):
    # raise NeedNewInternalsException(new_cart_coords)

    log(logger, "")

    # Return the difference between the new cartesian coordinates that yield
    # the desired internal coordinates and the old cartesian coordinates.
    cart_step = new_cart_coords - old_cart_coords
    return new_prim_ints, cart_step, backtransform_failed
