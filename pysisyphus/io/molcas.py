import h5py
import numpy as np

from pysisyphus.wavefunction.shells import Shell
from pysisyphus.wavefunction.shells_molcas import MolcasShells


def shells_from_molcas_h5(h5_fn, **kwargs):
    with h5py.File(h5_fn) as handle:
        # shape (nprims, 2); (exponent, contr. coeff) in every row
        primitives = handle["PRIMITIVES"][:]
        # shape (nprims, 3); (center_ind, L, shell_id) in every row
        primitive_ids = handle["PRIMITIVE_IDS"][:]
        # shape (ncenters, )
        atnums = handle["CENTER_ATNUMS"][:]
        # shape (ncenters, ); unused but checked to be the same as 'atnums'
        charges = handle["CENTER_CHARGES"][:]
        # shape (ncenters, 3); center coordinates in Bohr
        coords3d = handle["CENTER_COORDINATES"][:]

    # As I never thought about cases where these two arrays aren't the same
    # we check that they are the same! So we can deal with such cases when they
    # actually appear.
    np.testing.assert_allclose(atnums.astype(float), charges)

    shells = list()
    exps = list()
    coeffs = list()
    nprims = len(primitives)
    last_index = nprims - 1
    for i, prim_id in enumerate(primitive_ids):
        prim_id = tuple(prim_id)
        center_ind, L, _ = prim_id
        exponent, coeff = primitives[i]
        exps.append(exponent)
        coeffs.append(coeff)

        # Look ahead to see if we are at the end or if a new shell will begin.
        # In either case we create a new Shell object from the stored data and insert
        # it into the shells list. Then we reset the 'exps' and 'coeffs' list and starting
        # filling the next shell.
        if (i == last_index) or (prim_id != tuple(primitive_ids[i + 1])):
            ind0 = center_ind - 1
            shell = Shell(
                L,
                coords3d[ind0],
                coeffs,
                exps,
                ind0,
                atnums[ind0],
            )
            shells.append(shell)
            # Reset exponent & contraction coefficient lists
            exps = list()
            coeffs = list()
    # Outside of loop over primitives

    return MolcasShells(shells, **kwargs)
