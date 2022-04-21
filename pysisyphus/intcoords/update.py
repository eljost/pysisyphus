import numpy as np

from pysisyphus.config import BEND_MIN_DEG, LB_MIN_DEG
from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.eval import eval_primitives
from pysisyphus.intcoords.exceptions import NeedNewInternalsException
from pysisyphus.intcoords.valid import bend_valid, dihedral_valid
from pysisyphus.intcoords.PrimTypes import Bends, Dihedrals, Rotations


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
    bend_inds,
    check_dihedrals=False,
    check_bends=False,
    bend_min_deg=BEND_MIN_DEG,
    bend_max_deg=LB_MIN_DEG,
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

    invalid_inds = list()
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

    if check_bends and len(bend_inds) > 0:
        bends = [prim_internals[i] for i in bend_inds]
        are_valid = [
            bend_valid(new_coords3d, prim.inds, bend_min_deg, bend_max_deg)
            for prim in bends
        ]
        first_bend = bend_inds[0]
        invalid_inds = [
            i + first_bend for i, is_valid in enumerate(are_valid) if not is_valid
        ]

    if len(invalid_inds) > 0:
        invalid_prims = [primitives[i] for i in invalid_inds]
        invalid_msg = ", ".join([str(tp) for tp in invalid_prims])
        log(logger, "Internal coordinate(s) became invalid! Need new internal coordinates!")
        log(logger, f"Invalid primitives: {invalid_msg}")
        raise NeedNewInternalsException(
            new_coords3d, invalid_inds=invalid_inds, invalid_prims=invalid_prims
        )

    return prim_internals


def inds_from_prim_types(typed_prims, prim_types):
    inds = [
        i for i, (prim_type, *inds) in enumerate(typed_prims) if prim_type in prim_types
    ]
    return inds


def transform_int_step(
    int_step,
    old_cart_coords,
    cur_internals,
    Bt_inv_prim,
    primitives,
    typed_prims=None,
    dihedral_inds=None,
    rotation_inds=None,
    bend_inds=None,
    check_dihedrals=False,
    check_bends=False,
    bend_min_deg=BEND_MIN_DEG,
    bend_max_deg=LB_MIN_DEG,
    freeze_atoms=None,
    constrained_inds=None,
    update_constraints=False,
    cart_rms_thresh=1e-6,
    Bt_inv_prim_getter=None,
    max_cycles=25,
    logger=None,
):
    """Transformation is done in primitive internals, so int_step must be given
    in primitive internals and not in DLC!"""

    # If an iterable of typed prims is given we can derive bend/dihedral/rotation
    # indices from them.
    if typed_prims is not None:
        if dihedral_inds is None:
            dihedral_inds = inds_from_prim_types(typed_prims, Dihedrals)

        if bend_inds is None:
            bend_inds = inds_from_prim_types(typed_prims, Bends)

        if rotation_inds is None:
            rotation_inds = inds_from_prim_types(typed_prims, Rotations)

    if freeze_atoms is None:
        freeze_atoms = list()

    if constrained_inds is None:
        constrained_inds = list()

    freeze_atoms = np.array(freeze_atoms, dtype=int)
    new_cart_coords = old_cart_coords.copy()
    remaining_int_step = int_step
    target_internals = cur_internals + int_step
    # When we want to update the constraints we use the target primitive internals,
    if update_constraints:
        constrained_vals = target_internals[constrained_inds]
    # otherwise we try to stay with the original constrained values.
    else:
        constrained_vals = cur_internals[constrained_inds]

    def backtransform(remaining_int_step, Bt_inv_prim):
        """Separate function so it can be a called after the main loop
        finished to re-enforce the constraints."""
        nonlocal new_cart_coords

        cart_step = Bt_inv_prim.T.dot(remaining_int_step)
        # Remove step from frozen atoms.
        cart_step.reshape(-1, 3)[freeze_atoms] = 0.0
        cart_rms = np.sqrt(np.mean(cart_step**2))
        # Update cartesian coordinates
        new_cart_coords += cart_step
        # Determine new internal coordinates
        new_prim_ints = update_internals(
            new_cart_coords.reshape(-1, 3),
            old_internals,
            primitives,
            dihedral_inds,
            rotation_inds,
            bend_inds,
            check_dihedrals=check_dihedrals,
            check_bends=check_bends,
            bend_min_deg=bend_min_deg,
            bend_max_deg=bend_max_deg,
            logger=logger,
        )
        new_internals = [prim.val for prim in new_prim_ints]
        return new_prim_ints, new_internals, cart_rms

    last_rms = 9999
    old_internals = cur_internals
    backtransform_failed = True
    for i in range(max_cycles):
        if Bt_inv_prim_getter is not None:
            Bt_inv_prim = Bt_inv_prim_getter(new_cart_coords)
            log(logger, f"Recalculated (B^T)^+ in microcycle {i}")

        new_prim_ints, new_internals, cart_rms = backtransform(
            remaining_int_step, Bt_inv_prim
        )
        remaining_int_step = target_internals - new_internals
        internal_rms = np.sqrt(np.mean(remaining_int_step**2))
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

    if len(constrained_inds) > 0:
        for j in range(max_cycles):
            cur_constrained_vals = np.array(new_internals)[constrained_inds]
            diff = constrained_vals - cur_constrained_vals
            if any(np.abs(diff) <= 1e-5):
                break
            remaining_int_step = np.zeros_like(remaining_int_step)
            remaining_int_step[constrained_inds] = diff
            new_prim_ints, new_internals, _ = backtransform(
                remaining_int_step, Bt_inv_prim
            )
        if j > 0:
            log(logger, f"Re-enforced constraints in {j} additional cycle(s).")

    log(logger, "")

    # Return the difference between the new cartesian coordinates that yield
    # the desired internal coordinates and the old cartesian coordinates.
    cart_step = new_cart_coords - old_cart_coords
    return new_prim_ints, cart_step, backtransform_failed
