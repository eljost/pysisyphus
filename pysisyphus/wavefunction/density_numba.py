# [1] https://doi.org/10.1063/1.469408
#     Efficient molecular numerical integration schemes
#     Treutler, Ahlrichs, 1995


import math

import numba
import numpy as np


NORM2 = 0.5773502691896258
NORM3 = 0.2581988897471611
NORM4 = 0.09759000729485332
NORM22 = 0.3333333333333333


@numba.jit(
    nopython=True,
    cache=True,
    nogil=True,
)
def cart_gto3d_rel(La, axs, das, RA, RA2, result, cutoff=-40.0):
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
    shells, coords3d, P, precontr, rho, blk_size: int = 100, thresh: float = 1e-8
):
    """Density evaluation on a grid using numba.

    As (briefly) outlined in Section III C in [1]."""
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
        for i, R in enumerate(blk_coords3d):
            cur_center_ind = -1
            for shell in shells:
                La, A, center_ind, da, ax, a_ind, a_size = shell.as_tuple()
                # Only recompute distance to grid point when we are at a new center
                if center_ind != cur_center_ind:
                    RA = R - A
                    RA2 = RA**2
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
