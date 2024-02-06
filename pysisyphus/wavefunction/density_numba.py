# [1] https://doi.org/10.1063/1.469408
#     Efficient molecular numerical integration schemes
#     Treutler, Ahlrichs, 1995


import math

import numba
import numpy as np


CUTOFF = -36
NORM2 = 0.5773502691896258
NORM3 = 0.2581988897471611
NORM4 = 0.09759000729485332
NORM22 = 0.3333333333333333


@numba.jit(
    nopython=True,
    cache=True,
    nogil=True,
)
def cart_gto3d_rel(La, axs, das, RA, RA2, result, cutoff=CUTOFF):
    result[:] = 0.0
    dx, dy, dz = RA
    dx2, dy2, dz2 = RA2
    dist = dx2 + dy2 + dz2
    if La > 2:
        dx3 = dx2 * dx
        dy3 = dy2 * dy
        dz3 = dz2 * dz
    if La > 3:
        dx4 = dx2 * dx2
        dy4 = dy2 * dy2
        dz4 = dz2 * dz2

    for ax, da in zip(axs, das):
        exparg = -ax * dist
        if exparg < cutoff:
            continue
        expterm = da * math.exp(exparg)
        # s-Orbital
        if La == 0:
            result[0] += expterm
        # p-Orbital
        elif La == 1:
            result[0] += dx * expterm
            result[1] += dy * expterm
            result[2] += dz * expterm
        # d-Orbital
        elif La == 2:
            result[0] += NORM2 * dx2 * expterm
            result[1] += dx * dy * expterm
            result[2] += dx * dz * expterm
            result[3] += NORM2 * dy2 * expterm
            result[4] += dy * dz * expterm
            result[5] += NORM2 * dz2 * expterm
        # f-Orbital
        elif La == 3:
            result[0] += NORM3 * dx3 * expterm
            result[1] += NORM2 * dx2 * dy * expterm
            result[2] += NORM2 * dx2 * dz * expterm
            result[3] += NORM2 * dx * dy2 * expterm
            result[4] += dx * dy * dz * expterm
            result[5] += NORM2 * dx * dz2 * expterm
            result[6] += NORM3 * dy3 * expterm
            result[7] += NORM2 * dy2 * dz * expterm
            result[8] += NORM2 * dy * dz2 * expterm
            result[9] += NORM3 * dz3 * expterm
        # g-Orbital
        elif La == 4:
            result[0] += NORM4 * dx4 * expterm
            result[1] += NORM3 * dx3 * dy * expterm
            result[2] += NORM3 * dx3 * dz * expterm
            result[3] += NORM22 * dx2 * dy2 * expterm
            result[4] += NORM2 * dx2 * dy * dz * expterm
            result[5] += NORM22 * dx2 * dz2 * expterm
            result[6] += NORM3 * dx * dy3 * expterm
            result[7] += NORM2 * dx * dy2 * dz * expterm
            result[8] += NORM2 * dx * dy * dz2 * expterm
            result[9] += NORM3 * dx * dz3 * expterm
            result[10] += NORM4 * dy4 * expterm
            result[11] += NORM3 * dy3 * dz * expterm
            result[12] += NORM22 * dy2 * dz2 * expterm
            result[13] += NORM3 * dy * dz3 * expterm
            result[14] += NORM4 * dz4 * expterm
        else:
            # Could be changed to the explicit formula with a loop over a whole shell.
            result[:] = np.nan


@numba.jit(nopython=True, cache=True, parallel=True)
def eval_density(
    shells,
    coords3d,
    P,
    precontr,
    rho,
    blk_size: int = 100,
    thresh: float = 1e-8,
    accumulate: bool = False,
):
    """Density evaluation on a grid using numba.

    As (briefly) outlined in Section III C in [1]."""

    # Skip zeroing the density array when we accumulate.
    if not accumulate:
        rho[:] = 0.0
    # The grid will be split into blocks of a certain size
    nblks = int(np.ceil(len(coords3d) / blk_size))
    cart_size = sum([shell.cart_size() for shell in shells])
    sph_size = sum([shell.sph_size() for shell in shells])

    # Loop over all blocks
    for blk in numba.prange(nblks):
        blk_start = blk * blk_size
        blk_end = blk_start + blk_size
        blk_coords3d = coords3d[blk_start:blk_end]
        # As the number of grid points may not be a multiple of blk_size we have to
        # calculate the actual block size.
        cur_blk_size = len(blk_coords3d)

        chis_cart = np.zeros((cur_blk_size, cart_size))
        # The update to numba 0.59 made the declarations of RA and RA2 necessary.
        # The actual values of RA and RA2 will always be computed in the first
        # cycle of the loop over all shells, as 'cur_center_ind' will be NaN.
        RA = np.empty(3)
        RA2 = np.empty(3)
        for i, R in enumerate(blk_coords3d):
            cur_center_ind = np.nan
            for shell in shells:
                La, A, center_ind, da, ax, a_ind, a_size = shell.as_tuple()
                # Only recompute distance to grid point when we are at a new center
                if center_ind != cur_center_ind:
                    RA[:] = R - A
                    RA2[:] = RA**2
                    cur_center_ind = center_ind
                cart_gto3d_rel(
                    La, ax, da, RA, RA2, chis_cart[i, a_ind : a_ind + a_size]
                )

        # Convert to spherical basis functions and permute order
        chis = chis_cart @ precontr

        # Calculate mean values of basis functions in the current block
        chis_mean = np.empty(sph_size)
        for i in range(sph_size):
            chis_mean[i] = np.abs(chis[:, i]).mean()
        scratch = np.zeros((cur_blk_size, sph_size))

        # Keep track of which basis functions contributed
        nus_contribed = list()
        for nu in range(sph_size):
            for mu in range(nu, sph_size):
                factor = 1.0 if nu == mu else 2.0
                Pnm = P[nu, mu]
                contrib = factor * Pnm * chis_mean[nu] * chis_mean[mu]
                if abs(contrib) >= thresh:
                    scratch[:, nu] += factor * Pnm * chis[:, mu]
                    nus_contribed.append(nu)
        nus_contribed = sorted(set(nus_contribed))
        for nu in nus_contribed:
            rho[blk_start : blk_start + cur_blk_size] += scratch[:, nu] * chis[:, nu]


@numba.jit(
    nopython=True,
    cache=True,
    nogil=True,
)
def prefacts(La, RA, result):
    dx, dy, dz = RA
    if La > 1:
        dx2 = dx * dx
        dy2 = dy * dy
        dz2 = dz * dz
    if La > 2:
        dx3 = dx2 * dx
        dy3 = dy2 * dy
        dz3 = dz2 * dz
    if La > 3:
        dx4 = dx2 * dx2
        dy4 = dy2 * dy2
        dz4 = dz2 * dz2

    if La == 0:
        result[0] = 1.0
    # p-Orbital
    elif La == 1:
        result[0] = dx
        result[1] = dy
        result[2] = dz
    # d-Orbital
    elif La == 2:
        result[0] = NORM2 * dx2
        result[1] = dx * dy
        result[2] = dx * dz
        result[3] = NORM2 * dy2
        result[4] = dy * dz
        result[5] = NORM2 * dz2
    # f-Orbital
    elif La == 3:
        result[0] = NORM3 * dx3
        result[1] = NORM2 * dx2 * dy
        result[2] = NORM2 * dx2 * dz
        result[3] = NORM2 * dx * dy2
        result[4] = dx * dy * dz
        result[5] = NORM2 * dx * dz2
        result[6] = NORM3 * dy3
        result[7] = NORM2 * dy2 * dz
        result[8] = NORM2 * dy * dz2
        result[9] = NORM3 * dz3
    # g-Orbital
    elif La == 4:
        result[0] = NORM4 * dx4
        result[1] = NORM3 * dx3 * dy
        result[2] = NORM3 * dx3 * dz
        result[3] = NORM22 * dx2 * dy2
        result[4] = NORM2 * dx2 * dy * dz
        result[5] = NORM22 * dx2 * dz2
        result[6] = NORM3 * dx * dy3
        result[7] = NORM2 * dx * dy2 * dz
        result[8] = NORM2 * dx * dy * dz2
        result[9] = NORM3 * dx * dz3
        result[10] = NORM4 * dy4
        result[11] = NORM3 * dy3 * dz
        result[12] = NORM22 * dy2 * dz2
        result[13] = NORM3 * dy * dz3
        result[14] = NORM4 * dz4
    else:
        # Could be changed to the explicit formula with a loop over a whole shell.
        result[:] = np.nan


@numba.jit(nopython=True, cache=True)
def cart_size(L):
    return (L + 2) * (L + 1) // 2


@numba.jit(nopython=True, cache=True)
def eval_prim_density(Ls_inds, primdata, coords3d, P, switch, rho, cutoff=CUTOFF):
    # Initialize density array
    rho[:] = 0.0

    # Allocate prefactor arrays once using the maximum L values
    Lmax = np.max(Ls_inds[:, 0])
    prefacts_a = np.empty(cart_size(Lmax))
    prefacts_b = np.empty(cart_size(Lmax))

    nprims = len(primdata)
    # Loop over primitive pairs/first primitive
    for ai in range(nprims):
        da = primdata[ai, 0]
        ax = primdata[ai, 1]
        A = primdata[ai, 2:5]
        La, cart_index_a = Ls_inds[ai]
        cart_size_a = cart_size(La)

        # Loop over second primitive
        for bi in range(ai, nprims):
            db = primdata[bi, 0]
            bx = primdata[bi, 1]

            # Only carry out numerical integration for diffuse total exponents,
            # i.e, small total exponents px below the 'switch'-threshold.
            if ax + bx >= switch:
                continue

            B = primdata[bi, 2:5]
            Lb, cart_index_b = Ls_inds[bi]
            cart_size_b = cart_size(Lb)

            # Total exponent
            px = ax + bx
            # Reduced exponent
            mux = ax * bx / px
            # Skip primitive pair when exp-argument for pre-exponential factor
            # is too small.
            if -mux * np.sum((A - B) ** 2) <= cutoff:
                continue

            # Take symmetry into account, as we only loop over unique primitive
            # combinations.
            factor = 2.0
            if ai == bi:
                factor = 1.0

            # Pre-exponential factor w/ contraction coefficients and symmetry factor
            K = da * db * factor * np.exp(-mux * np.sum((A - B) ** 2))
            # Center-of-charge coordinate
            Povlp = (ax * A + bx * B) / px

            # Loop over grid points
            for i, R in enumerate(coords3d):
                RP = R - Povlp
                RP2 = RP**2
                # Check if exp-argument is below the threshold, if so we skip
                # this grid point.
                if -px * np.sum(RP2) <= cutoff:
                    continue
                Pexp = np.exp(-px * np.sum(RP2))

                # Determine angular momentum dependent prefactors of Gaussian
                # overlap distribution for both (primitive) shells.
                prefacts(La, R - A, prefacts_a)
                prefacts(Lb, R - B, prefacts_b)

                # Contract everything with the appropriate density matrix entries.
                for nu in range(cart_size_a):
                    for mu in range(cart_size_b):
                        rho[i] += (
                            K
                            * Pexp
                            * P[cart_index_a + nu, cart_index_b + mu]
                            * prefacts_a[nu]
                            * prefacts_b[mu]
                        )
