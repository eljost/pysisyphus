import itertools as it

import numpy as np

from pysisyphus.intcoords.eval import eval_primitives
from pysisyphus.intcoords import Torsion


def update_internals(new_coords3d, old_internals, primitives, dihedral_inds):
    prim_internals = eval_primitives(new_coords3d, primitives)
    new_internals = [prim_int.val for prim_int in prim_internals]
    internal_diffs = np.array(new_internals) - old_internals

    dihedrals = [prim_internals[i] for i in dihedral_inds]
    dihedral_num = len(dihedrals)
    dihedral_diffs = internal_diffs[-dihedral_num:]

    # Find differences that are shifted by 2*pi
    shifted_by_2pi = np.abs(np.abs(dihedral_diffs) - 2*np.pi) < np.pi/2
    new_dihedrals = np.array([dihed.val for dihed in dihedrals])
    if any(shifted_by_2pi):
        new_dihedrals[shifted_by_2pi] -= 2*np.pi * np.sign(dihedral_diffs[shifted_by_2pi])
        # Update values
        for dihed, new_val in zip(dihedrals, new_dihedrals):
            dihed.val = new_val

    return prim_internals


def transform_int_step(int_step, old_cart_coords, cur_internals, B_prim, primitives,
                       cart_rms_thresh=1e-6, logger=None):
    """Transformation is done in primitive internals, so int_step must be given
    in primitive internals and not in DLC!"""

    def log(msg):
        if logger is not None:
            logger.debug(msg)

    new_cart_coords = old_cart_coords.copy()
    remaining_int_step = int_step
    target_internals = cur_internals + int_step

    Bt_inv_prim = np.linalg.pinv(B_prim.dot(B_prim.T)).dot(B_prim)
    dihedral_inds = np.array([i for i, primitive in enumerate(primitives)
                 if isinstance(primitive, Torsion)
    ])

    last_rms = 9999
    old_internals = cur_internals
    backtransform_failed = True
    for i in range(25):
        cart_step = Bt_inv_prim.T.dot(remaining_int_step)
        cart_rms = np.sqrt(np.mean(cart_step**2))
        # Update cartesian coordinates
        new_cart_coords += cart_step
        # Determine new internal coordinates
        new_prim_ints = update_internals(
            new_cart_coords.reshape(-1, 3),
            old_internals,
            primitives,
            dihedral_inds
        )
        new_internals = [prim.val for prim in new_prim_ints]
        remaining_int_step = target_internals - new_internals
        internal_rms = np.sqrt(np.mean(remaining_int_step**2))
        log(f"Cycle {i}: rms(Δcart)={cart_rms:1.4e}, rms(Δint.) = {internal_rms:1.5e}")

        # This assumes the first cart_rms won't be > 9999 ;)
        if (cart_rms < last_rms):
            # Store results of the conversion cycle for laster use, if
            # the internal-cartesian-transformation goes bad.
            best_cycle = (new_cart_coords.copy(), new_internals.copy())
            best_cycle_ind = i
        elif i != 0:
            # If the conversion somehow fails we fallback to the best previous step.
            log(f"Backconversion failed! Falling back to step {best_cycle_ind}")
            new_cart_coords, new_internals = best_cycle
            break
        else:
            raise Exception("Internal-cartesian back-transformation already "
                            "failed in the first step. Aborting!"
            )
        old_internals = new_internals

        last_rms = cart_rms
        if cart_rms < cart_rms_thresh:
            log(f"Internal->Cartesian transformation converged in {i} cycle(s)!")
            backtransform_failed = False
            break

    # if self.check_dihedrals and (not self.dihedrals_are_valid(new_cart_coords)):
        # raise NeedNewInternalsException(new_cart_coords)

    log("")

    # Return the difference between the new cartesian coordinates that yield
    # the desired internal coordinates and the old cartesian coordinates.
    cart_step = new_cart_coords - old_cart_coords
    return new_prim_ints, cart_step, backtransform_failed
