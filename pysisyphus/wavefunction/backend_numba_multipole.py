import numba
import numpy as np


_FUNC_DATA = {
    #              Le,   components
    "int1e_ovlp": (0, 0),
    "int1e_r": (1, 3),
    "int1e_rr": (2, 9),
}


def get_func_data(key):
    # Le, components
    return _FUNC_DATA[key]


@numba.jit(nopython=True, cache=True)
def canonical_order(L: int) -> np.ndarray:
    inds = np.zeros(((L + 2) * (L + 1) // 2, 3), dtype=np.int64)
    i = 0
    for j in range(L + 1):
        l = L - j
        for n in range(j + 1):
            m = j - n
            inds[i] = (l, m, n)
            i += 1
    return inds


@numba.jit(nopython=True, cache=True)
def factorial2(n: int) -> int:
    """Double factorial for positive integer arguments and 0 and -1."""
    if n == -1:
        return 1
    offset = n % 2
    result = 1
    for i in range(2 + offset, n + 1, 2):
        result *= i
    return result


@numba.jit(nopython=True, cache=True)
def lmn_factors(l: int, m: int, n: int) -> float:
    return 1 / np.sqrt(
        factorial2(2 * l - 1) * factorial2(2 * m - 1) * factorial2(2 * n - 1)
    )


@numba.jit(nopython=True, cache=True)
def multipole1d_(i, j, e, px, pa, pb, pr, base):
    """1d-multipole integral."""
    # Base case
    if (i < 0) or (j < 0) or (e < 0) or ((i == 0) and (j == 0) and (e == 0)):
        # return np.sqrt(np.pi / px)
        return base
    # Decrement bra
    elif i > 0:
        return pa * multipole1d(i - 1, j, e, px, pa, pb, pr, base) + 1 / (2 * px) * (
            (i - 1) * multipole1d(i - 2, j, e, px, pa, pb, pr, base)
            + j * multipole1d(i - 1, j - 1, e, px, pa, pb, pr, base)
            + e * multipole1d(i - 1, j, e - 1, px, pa, pb, pr, base)
        )
    # Decrement ket
    elif j > 0:
        return pb * multipole1d(i, j - 1, e, px, pa, pb, pr, base) + 1 / (2 * px) * (
            i * multipole1d(i - 1, j - 1, e, px, pa, pb, pr, base)
            + (j - 1) * multipole1d(i, j - 2, e, px, pa, pb, pr, base)
            + e * multipole1d(i, j - 1, e - 1, px, pa, pb, pr, base)
        )
    # Decrement multipole order
    # e > 0
    else:
        return pr * multipole1d(i, j, e - 1, px, pa, pb, pr, base) + 1 / (2 * px) * (
            i * multipole1d(i - 1, j, e - 1, px, pa, pb, pr, base)
            + j * multipole1d(i, j - 1, e - 1, px, pa, pb, pr, base)
            + (e - 1) * multipole1d(i, j, e - 2, px, pa, pb, pr, base)
        )


@numba.jit(nopython=True, nogil=True, cache=True)
def multipole1d(i, j, e, px, pa, pb, pr, base):
    """1d-multipole integral."""

    def vrr(i, j, e):
        return multipole1d(i, j, e, px, pa, pb, pr, base)

    # Base case
    if (i < 0) or (j < 0) or (e < 0) or ((i == 0) and (j == 0) and (e == 0)):
        # return np.sqrt(np.pi / px)
        return base
    # Decrement bra
    elif i > 0:
        return pa * vrr(i - 1, j, e) + 1 / (2 * px) * (
            (i - 1) * vrr(i - 2, j, e)
            + j * vrr(i - 1, j - 1, e)
            + e * vrr(i - 1, j, e - 1)
        )
    # Decrement ket
    elif j > 0:
        return pb * vrr(i, j - 1, e) + 1 / (2 * px) * (
            i * vrr(i - 1, j - 1, e)
            + (j - 1) * vrr(i, j - 2, e)
            + e * vrr(i, j - 1, e - 1)
        )
    # Decrement multipole order
    # e > 0
    else:
        return pr * vrr(i, j, e - 1) + 1 / (2 * px) * (
            i * vrr(i - 1, j, e - 1)
            + j * vrr(i, j - 1, e - 1)
            + (e - 1) * vrr(i, j, e - 2)
        )


@numba.jit(nopython=True, cache=True, nogil=True, fastmath=True)
def multipole3d(La, Lb, Le, axs, das, A, bxs, dbs, B, R, exp_thresh=-36.0):
    """3d-multipole integral."""
    # Angular momenta of the different shells
    lmns_a = canonical_order(La)
    lmns_b = canonical_order(Lb)
    lmns_e = canonical_order(Le)
    na = len(lmns_a)
    nb = len(lmns_b)
    ne = len(lmns_e)
    # Final integrals
    integrals = np.zeros((na, nb, ne))

    # Construct angular momenta dependent normalization factors. Actually they
    # should be precalculated somehow.
    lmn_norms = np.zeros((na, nb))
    for i in range(na):
        la, ma, na_ = lmns_a[i]
        flmna = lmn_factors(la, ma, na_)
        for j in range(nb):
            lb, mb, nb_ = lmns_b[j]
            lmn_norms[i, j] = flmna * lmn_factors(lb, mb, nb_)

    # Number of primitives in bra and ket
    nprimsa = len(axs)
    nprimsb = len(bxs)
    AB = A - B
    AB2 = AB**2
    AB2sum = AB2.sum()

    # Determine most diffuse exponent pair in both shells and calulcate
    # the associated exp-argument. When this is already very small then
    # we skip the whole shell pair.
    ax_min = axs.min()
    bx_min = bxs.min()
    min_exp_arg = -(ax_min * bx_min) / (ax_min + bx_min) * AB2sum
    if min_exp_arg <= exp_thresh:
        return integrals

    # Loop over pairs of primitives
    for a in range(nprimsa):
        ax = axs[a]
        da = das[a]
        for b in range(nprimsb):
            bx = bxs[b]
            dadb = da * dbs[b]

            px = ax + bx
            mux = ax * bx / px
            exp_arg = -mux * AB2sum
            # Skip primitive pair when exp-argument is very small
            if exp_arg <= exp_thresh:
                continue
            K = np.exp(exp_arg)
            P = (ax * A + bx * B) / px
            PA = P - A
            PB = P - B
            PR = P - R
            base = np.sqrt(np.pi / px)
            # Loop over triples of angular momenta
            for i in range(na):
                lmna = lmns_a[i]
                for j in range(nb):
                    lmnb = lmns_b[j]
                    for k in range(ne):
                        lmne = lmns_e[k]
                        # Build up x-, y- and z-terms
                        tmp = 1.0
                        for m in range(3):
                            tmp *= multipole1d(
                                lmna[m], lmnb[m], lmne[m], px, PA[m], PB[m], PR[m], base
                            )
                        integrals[i, j, k] += dadb * K * tmp

    # Apply lmn-dependent basis function normalization
    for i in range(na):
        for j in range(nb):
            integrals[i, j] *= lmn_norms[i, j]
    return integrals


@numba.jit(parallel=True, nopython=True, cache=True)
def get_multipole_ints_cart_numba(
    Le,
    R,
    shells_a,
    shells_b,
    symmetric,
):
    components = 2 * Le + 1

    tot_size_a = 0
    for shell in shells_a:
        tot_size_a += shell.size

    tot_size_b = 0
    for shell in shells_b:
        tot_size_b += shell.size

    # Allocate final integral array
    integrals = np.zeros((tot_size_a, tot_size_b, components))
    shells_b = shells_a
    nshells_a = len(shells_a)
    nshells_b = len(shells_b)

    # Start parallel loop over contracted gaussians in shells_a
    for i in numba.prange(nshells_a):
        shell_a = shells_a[i]
        La, A, _, das, axs, indexa, sizea = shell_a.as_tuple()
        slicea = slice(indexa, indexa + sizea)

        # Start loop over contracted gaussians in shells_b
        if not symmetric:
            i = 0
        for j in range(i, nshells_b):
            shell_b = shells_b[j]
            Lb, B, _, dbs, bxs, indexb, sizeb = shell_b.as_tuple()
            sliceb = slice(indexb, indexb + sizeb)

            result = multipole3d(La, Lb, Le, axs, das, A, bxs, dbs, B, R)
            integrals[slicea, sliceb, :] = result

            if symmetric and (i != j):
                for k in range(indexa, indexa + sizea):
                    for l in range(indexb, indexb + sizeb):
                        integrals[l, k, :] = integrals[k, l, :]
        # End loop over contracted gaussians in shells_b
    # End loop over contracted gaussians in shells_a
    return integrals


def get_multipole_ints_cart(
    Le: int,
    R: np.ndarray,
    shellstructs_a,
    shellstructs_b=None,
):
    symmetric = shellstructs_b is None
    if symmetric:
        shellstructs_b = shellstructs_a

    integrals = get_multipole_ints_cart_numba(
        Le,
        R,
        shellstructs_a,
        shellstructs_b,
        symmetric,
    )
    if integrals.shape[2] == 1:
        integrals = np.squeeze(integrals, axis=2)
    return integrals


def get_1el_ints_cart(shells, func_dict, shells_b, **kwargs):
    R = kwargs.get("R", np.zeros(3))
    Le = func_dict
    return get_multipole_ints_cart(Le, R, shells, shells_b)


def get_multipole_ints_sph(Le: int, R: np.ndarray, shells_a, shells_b=None):
    """This function expects pysisyphus.Shells not Shellstructs"""
    shellstructs_a = shells_a.numba_shellstructs
    if shells_b is not None:
        shellstructs_b = shells_b.numba_shellstructs
    else:
        shellstructs_b = shellstructs_a

    integrals_cart = get_multipole_ints_cart(Le, R, shellstructs_a, shellstructs_b)

    c2s_coeffs_a = shells_a.reorder_c2s_coeffs
    if shells_b is not None:
        c2s_coeffs_b = shells_b.reorder_c2s_coeffs
    else:
        c2s_coeffs_b = c2s_coeffs_a

    int_matrix_sph = np.einsum(
        "ij,jk...,kl->il...",
        c2s_coeffs_a,
        integrals_cart,
        c2s_coeffs_b.T,
        optimize="greedy",
    )
    return int_matrix_sph
