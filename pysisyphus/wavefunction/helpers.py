from typing import Literal, Union

import numpy as np

from enum import IntEnum


L_MAP = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
}
L_SIZE = {l: (l + 1) * (l + 2) // 2 for l in L_MAP.values()}

Ls = Literal["s", "p", "d", "f", "g", "h"]
L_Inp = Union[int, Ls]


def get_l(l_inp: L_Inp) -> int:
    """Convert shell label to angular moment quantum number l."""
    try:
        l = L_MAP[l_inp.lower()]
    except (KeyError, AttributeError):
        l = int(l_inp)
    return l


def get_shell_shape(*Ls):
    return [L_SIZE[L] for L in Ls]


class BFType(IntEnum):
    CARTESIAN = 1
    PURE_SPHERICAL = 2


def canonical_order(L):
    inds = list()
    for i in range(L + 1):
        l = L - i
        for n in range(i + 1):
            m = i - n
            inds.append((l, m, n))
    return inds


def cca_order(l):
    """Same as canonical_order()."""
    inds = [
        (l, 0, 0),
    ]
    a = l
    b = c = 0
    for _ in range(((l + 1) * (l + 2) // 2) - 1):
        if c < l - a:
            b -= 1
            c += 1
        else:
            a -= 1
            c = 0
            b = l - a
        inds.append((a, b, c))
    return inds


def cart2ang_moms(cart):
    """Translate 'xxx' to (3, 0, 0) etc."""

    ang_moms = np.zeros(3, dtype=int)
    for char in cart:
        # x: 120, y: 121, z: 122
        ang_moms[ord(char.lower()) - 120] += 1
    return tuple(ang_moms)


def permut_for_order(cart_orders, pysis_order_func=cca_order):
    """Create permutation matrices to reorder Cartesian integrals."""
    Ps = dict()
    for cart_order in cart_orders:
        lengths = [len(_) for _ in cart_order]
        l0 = lengths[0]
        assert all([length == l0 for length in lengths])
        other_order = [cart2ang_moms(_) for _ in cart_order]
        pysis_order = pysis_order_func(l0)

        num = len(pysis_order)
        P = np.zeros((num, num), dtype=int)
        for i, po in enumerate(other_order):
            ind = pysis_order.index(po)
            P[i, ind] = 1
        Ps[l0] = P
    return Ps


def symmetric_orthogonalization(mat, S, thresh=1e-6):
    w, v = np.linalg.eigh(mat.T.dot(S).dot(mat))
    mask = w >= thresh
    w_isqrt = np.sqrt(w[mask])
    S_isqrt = (v[:, mask] / w_isqrt).dot(v[:, mask].T)
    return mat.dot(S_isqrt)
