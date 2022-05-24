import itertools as it
from typing import Tuple


import numpy as np
from numpy.typing import NDArray
import scipy as sp
from scipy.special import factorial2


from pysisyphus.config import L_MAX
from pysisyphus.elem_data import INV_ATOMIC_NUMBERS
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.helpers import canonical_order, get_l, get_shell_shape
from pysisyphus.wavefunction import dipoles3d, ovlps3d
from pysisyphus.wavefunction.cart2sph import cart2sph_coeffs


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
    def __init__(self, L, atomic_num, center, coeffs, exps, center_ind=None):
        self.L = get_l(L)
        self.atomic_num = atomic_num
        self.center = np.array(center, dtype=float)
        self.coeffs = np.array(coeffs, dtype=float)
        self.exps = np.array(exps, dtype=float)
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


def get_map(module, func_base_name):
    """Return dict that holds the different integrals functions."""
    return {
        (la, lb): getattr(module, f"{func_base_name}_{la}{lb}")
        for la, lb in it.product(Ls, Ls)
    }


Smap = get_map(ovlps3d, "ovlp3d")  # Overlap integrals
DPMmap = get_map(dipoles3d, "dipole3d")  # Dipole moments


def ovlp(la_tot, lb_tot, a, A, b, B):
    func = Smap[(la_tot, lb_tot)]
    return func(a, A, b, B)


def dpm(la_tot, lb_tot, a, A, b, B, C):
    func = DPMmap[(la_tot, lb_tot)]
    return func(a, A, b, B, C)


class Shells:
    sph_Ps = {l: np.eye(2 * l + 1) for l in range(L_MAX)}

    def __init__(self, shells, ordering="native"):
        self.shells = shells
        self.ordering = ordering
        assert ordering in ("pysis", "native")

    def __len__(self):
        return len(self.shells)

    def print_shells(self):
        for shell in self.shells:
            print(shell)

    @property
    def l_max(self):
        return max([shell.L for shell in self.shells])

    @property
    def atoms_coords3d(self):
        atoms = list()
        coords3d = list()
        center_inds = list()
        for shell in self.shells:
            center_ind = shell.center_ind
            if center_ind in center_inds:
                continue
            else:
                center_inds.append(center_ind)
            atom = INV_ATOMIC_NUMBERS[shell.atomic_num]
            atoms.append(atom)
            center = shell.center
            coords3d.append(center)
        coords3d = np.array(coords3d)
        return atoms, coords3d

    def from_basis(self, name, **kwargs):
        from pysisyphus.wavefunction.Basis import shells_with_basis

        atoms, coords3d = self.atoms_coords3d
        return shells_with_basis(atoms, coords3d, name=name, **kwargs)

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
        return self._ao_center_iter(lambda l: range(2 * l + 1))

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
                # Overlap over primitives, shape: (product(*shape), prims_a, prims_b)
                ss = ovlp(La, Lb, aa[:, None], A, bb[None, :], B)
                # Contraction coefficients
                # ss *= dA[None, :, None] * dB[None, None, :]
                ss *= dA[:, None] * dB[None, :]
                ss = ss.reshape(*shape, len(dA), len(dB))
                # Normalization
                ss *= normsa[:, None, :, None] * normsb[None, :, None, :]
                # Contract primitives, shape: shape
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

    def get_dipole_ints_cart(self, origin):
        shells_a = self
        shells_b = shells_a

        rows_x = list()
        rows_y = list()
        rows_z = list()
        for shell_a in shells_a.shells:
            La, A, dA, aa, normsa = shell_a.as_tuple()
            row_x = list()
            row_y = list()
            row_z = list()
            for shell_b in shells_b.shells:
                Lb, B, dB, bb, normsb = shell_b.as_tuple()
                shape = (*get_shell_shape(La, Lb), 3)
                dp = dpm(La, Lb, aa[:, None], A, bb[None, :], B, C=origin)
                dp *= dA[:, None] * dB[None, :]
                dp = dp.reshape(*shape, len(dA), len(dB))
                dp *= normsa[:, None, None, :, None] * normsb[None, :, None, None, :]
                # Shape: (*shape, 3)
                dp = dp.sum(axis=(3, 4)).reshape(shape)
                dp_x = dp[:, :, 0]
                dp_y = dp[:, :, 1]
                dp_z = dp[:, :, 2]
                row_x.append(dp_x)
                row_y.append(dp_y)
                row_z.append(dp_z)
            rows_x.append(row_x)
            rows_y.append(row_y)
            rows_z.append(row_z)
        DP_x = np.block(rows_x)
        DP_y = np.block(rows_y)
        DP_z = np.block(rows_z)
        return DP_x, DP_y, DP_z


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
