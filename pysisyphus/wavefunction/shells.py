from jinja2 import Template
import itertools as it
import textwrap
from typing import Tuple


import numpy as np
from numpy.typing import NDArray
import scipy as sp
from scipy.special import factorial2


from pysisyphus.config import L_MAX
from pysisyphus.elem_data import INV_ATOMIC_NUMBERS, nuc_charges_for_atoms
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.helpers import (
    canonical_order,
    get_l,
    get_shell_shape,
    permut_for_order,
)
from pysisyphus.wavefunction import coulomb3d, dipole3d, kinetic3d, ovlp3d

# from pysisyphus.wavefunction.devel_ints import coulomb3d, dipole3d, kinetic3d, ovlp3d
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
    def __init__(self, L, center, coeffs, exps, center_ind=None, atomic_num=None):
        self.L = get_l(L)
        self.center = np.array(center, dtype=float)
        self.coeffs = np.array(coeffs, dtype=float)
        self.exps = np.array(exps, dtype=float)
        assert self.coeffs.size == self.exps.size
        assert self.coeffs.shape == self.exps.shape
        self.center_ind = center_ind
        self.atomic_num = atomic_num

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

    def exps_coeffs_iter(self):
        return zip(self.exps, self.coeffs)

    @property
    def contr_depth(self):
        return self.coeffs.size

    def size(self):
        L = self.L
        return (L + 1) * (L + 2) // 2

    def __str__(self):
        try:
            center_str = f", at atom {self.center_ind}"
        except AttributeError:
            center_str = ""
        return f"Shell(L={self.L}, {self.contr_depth} pGTO{center_str})"

    def __repr__(self):
        return self.__str__()


Ls = list(range(L_MAX + 1))


def get_map(module, func_base_name, Ls_num=2):
    """Return dict that holds the different integrals functions."""
    func_map = dict()
    for ls in it.product(*[Ls for _ in range(Ls_num)]):
        ls_str = "".join([str(l) for l in ls])
        func_map[ls] = getattr(module, f"{func_base_name}_{ls_str}")
    return func_map


Smap = get_map(ovlp3d, "ovlp3d")  # Overlap integrals
Tmap = get_map(kinetic3d, "kinetic3d")  # Kinetic energy integrals
Vmap = get_map(coulomb3d, "coulomb3d")  # 1el Coulomb integrals
DPMmap = get_map(dipole3d, "dipole3d")  # Dipole moments integrals
# ERImap = get_map(eri, "eri", 4)  # Dipole moments integrals


def ovlp(la_tot, lb_tot, a, A, b, B):
    """Wrapper for overlap integrals."""
    func = Smap[(la_tot, lb_tot)]
    return func(a, A, b, B)


def kinetic(la_tot, lb_tot, a, A, b, B):
    """Wrapper for kinetic energy integrals."""
    func = Tmap[(la_tot, lb_tot)]
    return func(a, A, b, B)


def coulomb(la_tot, lb_tot, a, A, b, B, C):
    """Wrapper for 1el Coulomb integrals."""
    func = Vmap[(la_tot, lb_tot)]
    return func(a, A, b, B, C)


def dpm(la_tot, lb_tot, a, A, b, B, C):
    """Wrapper for linear moment integrals."""
    func = DPMmap[(la_tot, lb_tot)]
    return func(a, A, b, B, C)


# def eri(la_tot, lb_tot, lc_tot, ld_tot, a, A, b, B, c, C, d, D):
# """Wrapper for electron repulsion integrals."""
# func = ERImap[(la_tot, lb_tot, lc_tot, ld_tot)]
# return func(a, A, b, B, c, C, d, D)


class Shells:
    sph_Ps = {l: np.eye(2 * l + 1) for l in range(L_MAX)}

    def __init__(self, shells, ordering="native"):
        self.shells = shells
        self.ordering = ordering
        """
        'native' refers to the original ordering, as used in the QC program. The
        ordering will be reflected in the MO coefficient ordering. With 'native'
        the integrals calculated by pysisyphus must be reorderd, to match the native
        ordering of the MO coefficients.
        """
        assert ordering in ("pysis", "native")

        try:
            self.cart_Ps = permut_for_order(self.cart_order)
        except AttributeError:
            pass

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
            try:
                atom = INV_ATOMIC_NUMBERS[shell.atomic_num]
            # Use dummy atom when atomic_num is not set / None
            except KeyError:
                atom = "X"
            atoms.append(atom)
            center = shell.center
            coords3d.append(center)
        coords3d = np.array(coords3d)
        return atoms, coords3d

    @property
    def get_cart_bf_num(self):
        return sum([shell.size() for shell in self.shells])

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

    @staticmethod
    @file_or_str(".fchk")
    def from_fchk(text):
        # import here, to avoid cyclic imports
        from pysisyphus.io.fchk import shells_from_fchk

        shells = shells_from_fchk(text)
        return shells

    def center_shell_iter(self):
        sorted_shells = sorted(self.shells, key=lambda shell: shell.center_ind)
        return it.groupby(sorted_shells, key=lambda shell: shell.center_ind)

    def as_gau_gbs(self) -> str:
        def dfmt(num):
            return f"{num: 12.10e}".replace("e", "D")

        gbs_tpl = Template(
            """
        {% for center_ind, shells in center_shell_iter %}
            {{ center_ind+1 }} 0
        {%- for shell in shells %}
        {{ ("S", "P", "D", "F", "G")[shell.L] }}   {{ shell.exps.size}}    1.00 0.000000000000
            {%- for exp_, coeff in shell.exps_coeffs_iter() %}
            {{ dfmt(exp_) }} {{ dfmt(coeff) }}
            {%- endfor -%}
        {% endfor %}
         ****
        {%- endfor %}
        """
        )
        rendered = gbs_tpl.render(center_shell_iter=self.center_shell_iter(), dfmt=dfmt)
        rendered = textwrap.dedent(rendered)
        return rendered

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

    @property
    def P_cart(self):
        return sp.linalg.block_diag(*[self.cart_Ps[shell.L] for shell in self.shells])

    def get_1el_ints_cart(
        self, int_func, other=None, add_args=None, can_reorder: bool = True
    ) -> NDArray:
        shells_a = self
        shells_b = shells_a if other is None else other

        if add_args is None:
            add_args = {}

        rows = list()
        for shell_a in shells_a.shells:
            La, A, dA, aa, normsa = shell_a.as_tuple()
            row = list()
            for shell_b in shells_b.shells:
                Lb, B, dB, bb, normsb = shell_b.as_tuple()
                shape = get_shell_shape(La, Lb)
                # Integral over primitives, shape: (product(*shape), prims_a, prims_b)
                pints = int_func(La, Lb, aa[:, None], A, bb[None, :], B, **add_args)
                # Contraction coefficients
                pints *= dA[:, None] * dB[None, :]
                pints = pints.reshape(*shape, len(dA), len(dB))
                # Normalization
                pints *= normsa[:, None, :, None] * normsb[None, :, None, :]
                # Contract primitives, shape: shape
                cints = pints.sum(axis=(2, 3)).reshape(shape)
                row.append(cints)
            rows.append(row)
        int_matrix = np.block(rows)
        # Reordering will be disabled, when spherical integrals are desired. They
        # are reordered outside of this function. Reordering them already here
        # would mess up the results after the 2nd reordering.
        if can_reorder and self.ordering == "native":
            int_matrix = self.P_cart.dot(int_matrix).dot(shells_b.P_cart.T)
        return int_matrix

    def get_1el_ints_sph(
        self, int_func=None, int_matrix=None, other=None, **kwargs
    ) -> NDArray:
        if int_matrix is None:
            # Disallow reordering of the Cartesian integrals
            int_matrix = self.get_1el_ints_cart(
                int_func, other=other, can_reorder=False, **kwargs
            )

        shells_b = self if other is None else other
        C_a = self.cart2sph_coeffs
        C_b = shells_b.cart2sph_coeffs
        # Cartsian -> Spherical conversion
        int_matrix_sph = C_a.dot(int_matrix).dot(C_b.T)
        # Reorder to the native ordering of the respective program/format, e.g.
        # ORCA or .fchk.
        if self.ordering == "native":
            int_matrix_sph = self.P_sph.dot(int_matrix_sph).dot(shells_b.P_sph.T)
        return int_matrix_sph

    #####################
    # Overlap integrals #
    #####################

    def get_S_cart(self, other=None) -> NDArray:
        return self.get_1el_ints_cart(ovlp, other=other)

    def get_S_sph(self, other=None) -> NDArray:
        return self.get_1el_ints_sph(ovlp, other=other)

    @property
    def S_cart(self):
        return self.get_S_cart()

    @property
    def S_sph(self):
        return self.get_S_sph()

    ############################
    # Kinetic energy integrals #
    ############################

    def get_T_cart(self, other=None) -> NDArray:
        return self.get_1el_ints_cart(kinetic, other=other)

    def get_T_sph(self, other=None) -> NDArray:
        return self.get_1el_ints_sph(kinetic, other=other)

    @property
    def T_cart(self):
        return self.get_T_cart()

    @property
    def T_sph(self):
        return self.get_T_sph()

    ################################
    # Nuclear attraction integrals #
    ################################

    def get_V_cart(self, **kwargs) -> NDArray:
        atoms, coords3d = self.atoms_coords3d
        charges = nuc_charges_for_atoms(atoms)
        cart_bf_num = self.get_cart_bf_num
        V_nuc = np.zeros((cart_bf_num, cart_bf_num))
        # Loop over all cores and add the contributions
        for C, Z in zip(coords3d, charges):
            V_nuc += -Z * self.get_1el_ints_cart(
                coulomb,
                add_args={
                    "C": C,
                },
                **kwargs,
            )
        return V_nuc

    def get_V_sph(self) -> NDArray:
        V_cart = self.get_V_cart(can_reorder=False)
        return self.get_1el_ints_sph(int_func=None, int_matrix=V_cart)

    @property
    def V_cart(self):
        return self.get_V_cart()

    @property
    def V_sph(self):
        return self.get_V_sph()

    ###########################
    # Dipole moment integrals #
    ###########################

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
        return np.array((DP_x, DP_y, DP_z))

    def get_dipole_ints_sph(self, origin) -> NDArray:
        dp_cart = self.get_dipole_ints_cart(origin)
        shells_b = self
        C_a = self.cart2sph_coeffs
        C_b = shells_b.cart2sph_coeffs
        dp_sph = list()
        for dp in dp_cart:
            dp_ = C_a.dot(dp).dot(C_b.T)  # Cartsian -> Spherical conversion
            if self.ordering == "native":
                dp_ = self.P_sph.dot(dp_).dot(shells_b.P_sph.T)  # Reorder
            dp_sph.append(dp_)
        dp_sph = np.array(dp_sph)
        return dp_sph

    #######################################
    # Electron repulsion integrals - ERIs #
    #######################################

    def get_eri_cart(self):
        # shells_a = shells_b = shells_c = shells_d = self
        pass

    def __str__(self):
        return f"{self.__class__.__name__}({len(self.shells)} shells, ordering={self.ordering})"


class ORCAShells(Shells):
    sph_Ps = {
        0: [[1]],  # s
        1: [[0, 1, 0], [1, 0, 0], [0, 0, 1]],  # pz px py
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


class FCHKShells(Shells):
    cart_order = (
        ("",),
        ("x", "y", "z"),
        ("xx", "yy", "zz", "xy", "xz", "yz"),
        ("xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"),
        (
            "zzzz",
            "yzzz",
            "yyzz",
            "yyyz",
            "yyyy",
            "xzzz",
            "xyzz",
            "xyyz",
            "xyyy",
            "xxzz",
            "xxyz",
            "xxyy",
            "xxxz",
            "xxxy",
            "xxxx",
        ),
    )

    sph_Ps = {
        0: [[1]],  # s
        1: [[1, 0, 0], [0, 0, 1], [0, 1, 0]],  # px pz py
        2: [
            [0, 0, 1, 0, 0],  # dz²
            [0, 1, 0, 0, 0],  # dxz
            [0, 0, 0, 1, 0],  # dyz
            [1, 0, 0, 0, 0],  # dx² - y²
            [0, 0, 0, 0, 1],  # dxy
        ],
        # Similar to ORCA, but w/o the sign flips ;)
        3: [
            [0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0],
            [1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1],
        ],
        4: [
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
    }


class MoldenShells(Shells):
    sph_Ps = {
        0: [[1]],  # s
        1: [[1, 0, 0], [0, 0, 1], [0, 1, 0]],  # px, py, pz
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
