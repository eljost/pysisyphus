import warnings

import numpy as np


def bf_ind_key(bf_ind):
    """Key function to sort OpenMolcas basis functions."""
    center_ind, shell_ind, L, m = bf_ind
    if L == 1:
        key = (center_ind, L, m, shell_ind)
    else:
        key = (center_ind, L, -m, shell_ind)
    return key


def get_molcas_P_sph(shells, nbfs):
    warnings.warn(
        "'get_molcas_P_sph() assumes that all basis functions but the p-functions "
        "are spherical! p-functions are assumed to be Cartesian."
    )
    bf_inds = list()
    prev_key = None
    shell_ind = 0  # Define shell_ind to satisfy the linter
    for shell in shells:
        center_ind = shell.center_ind
        L = shell.L
        key = (center_ind, L)
        # Still at the same center and have the same angular momentum as before.
        # Increment the shell index by one.
        if prev_key and key == prev_key:
            shell_ind += 1
        # Indicate new shell by incremeting shell index
        else:
            shell_ind = 0
        for m in range(-L, L + 1):
            bf_ind = (center_ind, shell_ind, L, m)
            bf_inds.append(bf_ind)
        prev_key = key

    # After adding 1 (one) to the first two columns, the array
    # 'bf_inds_sorted' would correspond to BASIS_FUNCTION_IDS in the
    # Molcas HD5 files.
    # bf_inds_sorted = sorted(bf_inds, key=bf_ind_key)

    inds = sorted(range(len(bf_inds)), key=lambda i: bf_ind_key(bf_inds[i]))
    P_sph = np.zeros((nbfs, nbfs))
    P_sph[np.arange(nbfs), inds] = 1
    return P_sph
