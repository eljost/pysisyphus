from typing import Literal, Union, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.special import factorial2


from pysisyphus.helpers_pure import file_or_str


Ls = Literal["s", "p", "d", "f", "g", "h"]
L_Inp = Union[int, Ls]

L_MAP = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
}
L_SIZE = {l: (l + 1) * (l + 2) // 2 for l in L_MAP.values()}


def get_l(l_inp: L_Inp) -> int:
    """Convert shell label to angular moment quantum number l."""
    try:
        l = L_MAP[l_inp.lower()]
    except (KeyError, AttributeError):
        l = int(l_inp)
    return l


def get_shell_shape(La, Lb):
    return (L_SIZE[La], L_SIZE[Lb])


def canonical_order(L):
    inds = list()
    for i in range(L + 1):
        l = L - i
        for n in range(i + 1):
            m = i - n
            inds.append((l, m, n))
    return inds


def normalize(lmn: Tuple[int, int, int], coeffs: NDArray, exps: NDArray):
    """Norm of contracted GTO and norms of primitive GTOs.

    Based on a function, originally published by Joshua Goings here:

        https://joshuagoings.com/2017/04/28/integrals/
    """
    l, m, n = lmn
    L = l + m + n
    # self.norm is a list of length equal to number primitives
    # normalize primitives first (PGBFs)
    norm = np.sqrt(
        np.power(2, 2 * L + 1.5)
        * np.power(exps, L + 1.5)
        / factorial2(2 * l - 1)
        / factorial2(2 * m - 1)
        / factorial2(2 * n - 1)
        / np.power(np.pi, 1.5)
    )

    # now normalize the contracted basis functions (CGBFs)
    # Eq. 1.44 of Valeev integral whitepaper
    prefactor = (
        np.power(np.pi, 1.5)
        * factorial2(2 * l - 1)
        * factorial2(2 * m - 1)
        * factorial2(2 * n - 1)
        / np.power(2.0, L)
    )

    N = 0.0
    num_exps = len(exps)
    for ia in range(num_exps):
        for ib in range(num_exps):
            N += (
                norm[ia]
                * norm[ib]
                * coeffs[ia]
                * coeffs[ib]
                / np.power(exps[ia] + exps[ib], L + 1.5)
            )

    N *= prefactor
    N = np.power(N, -0.5)
    return N, norm


class Shell:
    def __init__(self, L, center, coeffs, exps, center_ind=None):
        self.L = get_l(L)
        self.center = np.array(center)
        self.coeffs = np.array(coeffs)
        self.exps = np.array(exps)
        assert self.coeffs.size == self.exps.size
        assert self.coeffs.shape == self.exps.shape
        self.center_ind = center_ind

        self._norms = None

    @property
    def norms(self):
        if self._norms is None:
            lmns = canonical_order(self.L)
            Ns = list()
            norms = list()
            for lmn in lmns:
                N, norm = normalize(lmn, self.coeffs, self.exps)
                Ns.append(N)
                norms.append(norm)
            self._norms = np.array(norms)
        return self._norms

    def as_tuple(self):
        return (self.L, self.center, self.coeffs, self.exps, self.norms)

    @property
    def contr_depth(self):
        return self.coeffs.size

    def size(self):
        L = self.L
        return (L + 1) * (L + 2) // 2

    def __str__(self):
        center_str = f", at atom {self.center_ind}" if self.center_ind else ""
        return f"Shell(L={self.L}, {self.contr_depth} pGTO(s){center_str})"

    def __repr__(self):
        return self.__str__()


import itertools as it
from pysisyphus.integrals import ovlps3d

Ls = (0, 1, 2)
Smap = {(la, lb): getattr(ovlps3d, f"ovlp3d_{la}{lb}") for la, lb in it.product(Ls, Ls)}


def ovlp(la_tot, lb_tot, a, A, b, B):
    func = Smap[(la_tot, lb_tot)]
    return func(a, A, b, B)


class Shells:
    def __init__(self, shells):
        self.shells = shells

    def __len__(self):
        return len(self.shells)

    @staticmethod
    @file_or_str(".in")
    def from_aomix(text):
        # import here, to avoid cyclic imports
        from pysisyphus.io.aomix import parse_aomix

        shells = parse_aomix(text)
        return shells

    def get_S_cart(self, other=None) -> NDArray:
        shells_a = self.shells
        shells_b = shells_a if other is None else other

        rows = list()
        for shell_a in shells_a:
            La, A, dA, aa, normsa = shell_a.as_tuple()
            row = list()
            for shell_b in shells_b:
                Lb, B, dB, bb, normsb = shell_b.as_tuple()
                shape = get_shell_shape(La, Lb)
                ss = ovlp(La, Lb, aa[:, None], A, bb[None, :], B)
                ss *= dA[None, :, None] * dB[None, None, :]
                ss = ss.reshape(*shape, len(dA), len(dB))
                ss *= normsa[:, None, :, None] * normsb[None, :, None, :]
                ss = ss.sum(axis=(2, 3)).reshape(shape)
                row.append(ss)
            rows.append(row)
        S = np.block(rows)
        return S

    @property
    def S_cart(self):
        return self.get_S_cart()
