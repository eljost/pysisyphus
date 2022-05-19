import itertools as it
from typing import Literal, Union, Tuple


import numpy as np
from numpy.typing import NDArray
import scipy as sp
from scipy.special import factorial2


from pysisyphus.config import L_MAX
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.helpers import canonical_order
from pysisyphus.wavefunction import ovlps3d
from pysisyphus.wavefunction.cart2sph import cart2sph_coeffs


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


def get_shell_shape(La, Lb):
    return (L_SIZE[La], L_SIZE[Lb])


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


Ls = list(range(L_MAX + 1))
Smap = {(la, lb): getattr(ovlps3d, f"ovlp3d_{la}{lb}") for la, lb in it.product(Ls, Ls)}


def ovlp(la_tot, lb_tot, a, A, b, B):
    func = Smap[(la_tot, lb_tot)]
    return func(a, A, b, B)


class Shells:
    def __init__(self, shells, ordering="native"):
        self.shells = shells
        self.ordering = ordering
        assert ordering in ("pysis", "native")

    def __len__(self):
        return len(self.shells)

    @property
    def l_max(self):
        return max([shell.L for shell in self.shells])

    @staticmethod
    @file_or_str(".in")
    def from_aomix(text):
        # import here, to avoid cyclic imports
        from pysisyphus.io.aomix import parse_aomix

        shells = parse_aomix(text)
        return shells

    @staticmethod
    @file_or_str(".json")
    def from_orca_json(text):
        # import here, to avoid cyclic imports
        from pysisyphus.io.orca import shells_from_json

        shells = shells_from_json(text)
        return shells
    
    def _ao_center_iter(self, func):
        i = 0
        for shell in self.shells:
            for _ in func(shell.L):
                yield shell.center_ind
                i += 0

    @property
    def cart_ao_centers(self):
        return self._ao_center_iter(canonical_order)

    @property
    def sph_ao_centers(self):
        return self._ao_center_iter(lambda l: range(2*l + 1))

    @property
    def cart2sph_coeffs(self):
        cart2sph = cart2sph_coeffs(self.l_max)
        C = sp.linalg.block_diag(*[cart2sph[shell.L] for shell in self.shells])
        return C

    @property
    def P_sph(self):
        return sp.linalg.block_diag(*[self.sph_Ps[shell.L] for shell in self.shells])

    def get_S_cart(self, other=None) -> NDArray:
        shells_a = self
        shells_b = shells_a if other is None else other

        rows = list()
        for shell_a in shells_a.shells:
            La, A, dA, aa, normsa = shell_a.as_tuple()
            row = list()
            for shell_b in shells_b.shells:
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

    def get_S_sph(self, other=None) -> NDArray:
        S_cart = self.get_S_cart(other=other)
        shells_b = self if other is None else other
        C_a = self.cart2sph_coeffs
        C_b = shells_b.cart2sph_coeffs
        S_sph = C_a.dot(S_cart).dot(C_b.T)  # Cartsian -> Spherical conversion
        if self.ordering == "native":
            S_sph = self.P_sph.dot(S_sph).dot(shells_b.P_sph.T)  # Reorder
        return S_sph

    @property
    def S_cart(self):
        return self.get_S_cart()

    @property
    def S_sph(self):
        return self.get_S_sph()


class ORCAShells(Shells):
    sph_Ps = {
        0: [[1]],  # s
        1: [[0, 1, 0], [1, 0, 0], [0, 0, 1]],  # p  pz px py
        2: [
            [0, 0, 1, 0, 0],  # dz²
            [0, 1, 0, 0, 0],  # dxz
            [0, 0, 0, 1, 0],  # dyz
            [1, 0, 0, 0, 0],  # dx² - y²
            [0, 0, 0, 0, 1],  # dxy
        ],
        3: [
            [0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0],
            [-1, 0, 0, 0, 0, 0, 0],  # ORCA, why you do this to me?
            [0, 0, 0, 0, 0, 0, -1],  # sign flip, as line abo ve
        ],
        4: [
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, -1, 0, 0, 0, 0, 0, 0, 0],  # sign flip
            [0, 0, 0, 0, 0, 0, 0, -1, 0],  # sign flip
            [-1, 0, 0, 0, 0, 0, 0, 0, 0],  # sign flip
            [0, 0, 0, 0, 0, 0, 0, 0, -1],  # sign flip
        ],
    }
