import numpy as np


from pysisyphus.wavefunction.ints import (
    coulomb3d,
    diag_quadrupole3d,
    dipole3d,
    kinetic3d,
    ovlp3d,
    quadrupole3d,
    int2c2e3d,
    int3c2e3d_sph,
)

_FUNC_DATA = {
    #              func_dict,     components
    "int1e_ovlp": (ovlp3d.ovlp3d, 0),
    "int1e_r": (dipole3d.dipole3d, 3),
    "int1e_rr": (quadrupole3d.quadrupole3d, 6),
    "int1e_drr": (diag_quadrupole3d.diag_quadrupole3d, 3),
    "int1e_kin": (kinetic3d.kinetic3d, 0),
    "int1e_nuc": (coulomb3d.coulomb3d, 0),
    "int2c2e": (int2c2e3d.int2c2e3d, 0),
    "int3c2e_sph": (int3c2e3d_sph.int3c2e3d_sph, 3),
}


def get_func_data(key):
    # func_dict, components
    return _FUNC_DATA[key]


def get_1el_ints_cart(
    shells_a,
    func_dict,
    shells_b=None,
    can_reorder=True,
    ordering="native",
    components=0,
    screen=False,
    screen_func=None,
    R=np.zeros(3),
    **kwargs,
):
    symmetric = shells_b is None
    shells_b = shells_a if symmetric else shells_b
    cart_bf_num_a = shells_a.cart_bf_num
    cart_bf_num_b = shells_b.cart_bf_num

    # components = 0 indicates, that a plain 2d matrix is desired.
    if is_2d := (components == 0):
        components = 1

    # Preallocate empty matrices and directly assign the calculated values
    integrals = np.zeros((components, cart_bf_num_a, cart_bf_num_b))

    # Dummy screen function that always returns True
    if (not screen) or (screen_func is None):
        screen_func = lambda *args: True

    for i, shell_a in enumerate(shells_a):
        La, A, da, ax, a_ind, a_size = shell_a.as_tuple()
        a_slice = slice(a_ind, a_ind + a_size)
        # When we don't deal with symmetric matrices we have to iterate over
        # all other basis functions in shells_b, so we reset i to 0.
        if not symmetric:
            i = 0
        a_min_exp = ax.min()  # Minimum exponent used for screening
        for shell_b in shells_b[i:]:
            Lb, B, db, bx, b_ind, b_size = shell_b.as_tuple()
            b_slice = slice(b_ind, b_ind + b_size)
            # Screen shell pair
            b_min_exp = bx.min()  # Minimum exponent used for screening
            R_ab = np.linalg.norm(A - B)  # Shell center distance
            if not screen_func(a_min_exp, b_min_exp, R_ab):
                continue
            integrals[:, a_slice, b_slice] = func_dict[(La, Lb)](
                ax[:, None],
                da[:, None],
                A,
                bx[None, :],
                db[None, :],
                B,
                # R=R,
                **kwargs,
            )
            # Fill up the lower tringular part of the matrix in the symmetric
            # case.
            # TODO: this may be (?) better done outside of this loop, after
            # everything is calculated...
            if symmetric:
                integrals[:, b_slice, a_slice] = np.transpose(
                    integrals[:, a_slice, b_slice], axes=(0, 2, 1)
                )

    # Return plain 2d array if components is set to 0, i.e., remove first axis.
    if is_2d:
        integrals = np.squeeze(integrals, axis=0)

    # Reordering will be disabled, when spherical integrals are desired. They
    # are reordered outside of this function. Reordering them already here
    # would mess up the results after the 2nd reordering.
    if can_reorder and ordering == "native":
        integrals = np.einsum(
            "ij,...jk,kl->...il",
            shells_a.P_cart,
            integrals,
            shells_b.P_cart.T,
            optimize="greedy",
        )
    return integrals


def get_3c2el_ints_cart(shells_a, shells_aux):
    """Cartesian 3-center-2-electron integrals.

    DO NOT USE THESE INTEGRALS AS THEY ARE RETURNED FROM THIS METHOD.
    These integrals utilize recurrence relations that are only valid,
    when the resulting Cartesian integrals are transformed into spherical
    integrals.

    Contrary to the general function 'get_1el_ints_cart', that supports
    different 'func_dict' arguments and cross-integrals between two
    different shells this function is less general. This function is
    restricted to '_3center2el_sph' and always uses only 1 set of shells
    and 1 set of auxiliary shells.
    """
    cart_bf_num_a = shells_a.cart_bf_num
    cart_bf_num_aux = shells_aux.cart_bf_num

    integrals = np.zeros((cart_bf_num_a, cart_bf_num_a, cart_bf_num_aux))
    func_dict = int3c2e3d_sph.int3c2e3d_sph

    for i, shell_a in enumerate(shells_a):
        La, A, da, ax, a_ind, a_size = shell_a.as_tuple()
        a_slice = slice(a_ind, a_ind + a_size)
        # As noted in the docstring, we iterate over pairs of shells_a
        for shell_b in shells_a[i:]:
            Lb, B, db, bx, b_ind, b_size = shell_b.as_tuple()
            b_slice = slice(b_ind, b_ind + b_size)
            for shell_c in shells_aux:
                Lc, C, dc, cx, c_ind, c_size = shell_c.as_tuple()
                c_slice = slice(c_ind, c_ind + c_size)
                integrals[a_slice, b_slice, c_slice] = func_dict[(La, Lb, Lc)](
                    ax[:, None, None],
                    da[:, None, None],
                    A,
                    bx[None, :, None],
                    db[None, :, None],
                    B,
                    cx[None, None, :],
                    dc[None, None, :],
                    C,
                )
                integrals[b_slice, a_slice, c_slice] = np.transpose(
                    integrals[a_slice, b_slice, c_slice], axes=(1, 0, 2)
                )
    return integrals
