import itertools as it
from math import sqrt, log, pi
from pathlib import Path
import textwrap
from typing import List, Literal, Tuple


from jinja2 import Template
from joblib import Memory
import numpy as np
from numpy.typing import NDArray
import scipy as sp
from scipy.special import factorial2


from pysisyphus.config import L_MAX, L_AUX_MAX
from pysisyphus.elem_data import (
    ATOMIC_NUMBERS,
    INV_ATOMIC_NUMBERS,
    nuc_charges_for_atoms,
)
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction.helpers import (
    canonical_order,
    get_l,
    get_shell_shape,
    permut_for_order,
)

from pysisyphus.wavefunction.ints import (
    coulomb3d,
    diag_quadrupole3d,
    dipole3d,
    gto3d,
    kinetic3d,
    ovlp3d,
    quadrupole3d,
    _2center2el3d,
    _3center2el3d,
    _3center2el3d_sph,
)

from pysisyphus.wavefunction.cart2sph import cart2sph_coeffs


def normalize(lmn: Tuple[int, int, int], coeffs: NDArray, exps: NDArray):
    """Norm of primitive GTOs and corresponding coefficients to normalize them.

    Based on a function, originally published by Joshua Goings here:

        https://joshuagoings.com/2017/04/28/integrals/
    """
    coeffs = np.atleast_2d(coeffs)
    L = sum(lmn)
    fact2l, fact2m, fact2n = [factorial2(2 * _ - 1) for _ in lmn]

    # Normalize primitives.
    norm = np.sqrt(
        np.power(2, 2 * L + 1.5)
        * np.power(exps, L + 1.5)
        / fact2l
        / fact2m
        / fact2n
        / np.power(np.pi, 1.5)
    )

    # Normalize the contracted basis functions (CGBFs)
    # Eq. 1.44 of Valeev integral whitepaper.
    prefactor = np.power(np.pi, 1.5) * fact2l * fact2m * fact2n / np.power(2.0, L)

    N = (coeffs.T * coeffs / np.power(exps[:, None] + exps[None, :], L + 1.5)).sum()
    N *= prefactor
    N = 1 / np.sqrt(N)
    return N, norm


def eval_pgtos(xyz, center, exponents, cart_powers):
    """Evaluate primititve Cartesian GTOs at points xyz."""

    xa, ya, za = (xyz - center).T
    # Indepdendent of Cartesian powers, but dependent on contraction degree
    exp_term = np.exp(-exponents[:, None] * (xa**2 + ya**2 + za**2))
    # Indepdendent of contraction degree, but dependent on Cartesian powers
    xpow, ypow, zpow = cart_powers.T
    prefac = xa ** xpow[:, None] * ya ** ypow[:, None] * za ** zpow[:, None]
    return prefac[:, None, :] * exp_term


def eval_shell(xyz, center, cart_powers, contr_coeffs, exponents):
    # pgtos shape is (nbf, npgto, nxyz)
    # that is
    # (number of basis funcs in shell, number of primitives, number of points)
    pgtos = eval_pgtos(xyz, center, exponents, cart_powers)
    # shape of contr_coeffs is (number of bfs in shell, number of primitives)
    pgtos *= contr_coeffs[..., None]
    return pgtos.sum(axis=1)


class Shell:
    def __init__(self, L, center, coeffs, exps, center_ind=None, atomic_num=None):
        self.L = get_l(L)
        self.center = np.array(center, dtype=float)  # (x, y, z), 1d array
        # Contraction coefficients, 2d array, shape (bfs in shell, number of primitives)
        self.org_coeffs = np.array(coeffs, dtype=float)
        self.coeffs = np.atleast_2d(self.org_coeffs.copy())
        self.exps = np.array(exps, dtype=float)  # Orbital exponents, 1d array
        assert self.coeffs.shape[1] == self.exps.size
        self.center_ind = int(center_ind)
        self.atomic_num = int(atomic_num)

        # Store original copy of contraction coefficients
        self.coeffs_org = self.coeffs.copy()
        # Compute norms and multiply them onto the contraction coefficients
        _, norm_coeffs = self.get_norms()
        self.coeffs = self.coeffs * norm_coeffs

    def get_norms(self, coeffs=None):
        """Shape (nbfs, nprimitives).

        For s- and p-orbitals all rows will be the same, but they will start
        to differ from d-orbitals on."""
        if coeffs is None:
            coeffs = self.coeffs
        lmns = canonical_order(self.L)
        Ns, norm_coeffs = zip(*[normalize(lmn, coeffs, self.exps) for lmn in lmns])
        return Ns, np.array(norm_coeffs)

    def as_tuple(self):
        return self.L, self.center, self.coeffs, self.exps

    def exps_coeffs_iter(self):
        return zip(self.exps, self.coeffs)

    @property
    def contr_depth(self):
        return self.coeffs.size

    @property
    def cart_powers(self):
        return np.array(canonical_order(self.L), dtype=int)

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


def get_map(module, func_base_name, Ls=(L_MAX, L_MAX)):
    """Return dict that holds the different integrals functions."""
    func_map = dict()
    L_ranges = [range(L + 1) for L in Ls]
    for ls in it.product(*L_ranges):
        ls_str = "".join([str(l) for l in ls])
        try:
            func_map[ls] = getattr(module, f"{func_base_name}_{ls_str}")
        except AttributeError:
            pass
    return func_map


CGTOmap = get_map(gto3d, "cart_gto3d", Ls=(L_MAX,))  # Cartesian GTO shells
Smap = get_map(ovlp3d, "ovlp3d")  # Overlap integrals
Tmap = get_map(kinetic3d, "kinetic3d")  # Kinetic energy integrals
Vmap = get_map(coulomb3d, "coulomb3d")  # 1el Coulomb integrals
DPMmap = get_map(dipole3d, "dipole3d")  # Dipole moments integrals
QPMmap = get_map(quadrupole3d, "quadrupole3d")  # Quadrupole moments integrals
DQPMmap = get_map(
    diag_quadrupole3d, "diag_quadrupole3d"
)  # Diagonal quadrupole moments integrals
_2c2elMap = get_map(_2center2el3d, "_2center2el3d", Ls=(L_AUX_MAX, L_AUX_MAX))
_3c2elMap = get_map(_3center2el3d, "_3center2el3d", Ls=(L_MAX, L_MAX, L_AUX_MAX))
_3c2elSphMap = get_map(
    _3center2el3d_sph, "_3center2el3d_sph", Ls=(L_MAX, L_MAX, L_AUX_MAX)
)


def cart_gto(l_tot, a, Xa, Ya, Za):
    """Wrapper for evaluation of cartesian GTO shells."""
    func = CGTOmap[(l_tot,)]
    return func(a, Xa, Ya, Za)


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


def dipole(la_tot, lb_tot, a, A, b, B, C):
    """Wrapper for linear moment integrals (dipole moment)."""
    func = DPMmap[(la_tot, lb_tot)]
    return func(a, A, b, B, C)


def quadrupole(la_tot, lb_tot, a, A, b, B, C):
    """Wrapper for quadratic moment integrals (quadrupole moment)."""
    func = QPMmap[(la_tot, lb_tot)]
    return func(a, A, b, B, C)


def diag_quadrupole(la_tot, lb_tot, a, A, b, B, C):
    """Wrapper for diagonal entries of quadratic moment integrals."""
    func = DQPMmap[(la_tot, lb_tot)]
    return func(a, A, b, B, C)


def _2center2electron(la_tot, lb_tot, a, A, b, B):
    """Wrapper for 2-center-2-electron integrals."""
    func = _2c2elMap[(la_tot, lb_tot)]
    return func(a, A, b, B)


def _3center2electron(la_tot, lb_tot, lc_tot, a, A, b, B, c, C):
    """Wrapper for 3-center-2-electron integrals."""
    func = _3c2elMap[(la_tot, lb_tot, lc_tot)]
    return func(a, A, b, B, c, C)


def _3center2electron_sph(la_tot, lb_tot, lc_tot, a, A, b, B, c, C):
    """Wrapper for 3-center-2-electron integrals."""
    func = _3c2elSphMap[(la_tot, lb_tot, lc_tot)]
    return func(a, A, b, B, c, C)


Ordering = Literal["native", "pysis"]


class Shells:
    sph_Ps = {l: np.eye(2 * l + 1) for l in range(L_MAX + 1)}

    def __init__(
        self,
        shells: List[Shell],
        screen: bool = True,
        cache: bool = True,
        cache_path: str = "./cache",
        ordering: Ordering = "native",
    ):
        self.shells = shells
        self.ordering = ordering
        self.screen = screen
        self.cache = cache
        self.cache_path = Path(cache_path)
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

        self.atoms, self.coords3d = self.atoms_coords3d
        # Precontract & store coefficients for reordering spherical basis functions
        # and converting them from Cartesian basis functions.
        self.reorder_c2s_coeffs = self.P_sph @ self.cart2sph_coeffs

        # Enable disk cache for 1el-integrals
        if self.cache:
            self.memory = Memory(self.cache_path, verbose=0)
            self.get_1el_ints_cart = self.memory.cache(self.get_1el_ints_cart)
            self.get_1el_ints_sph = self.memory.cache(self.get_1el_ints_sph)

    def __len__(self):
        return len(self.shells)

    def __getitem__(self, key):
        return self.shells[key]

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
            # Skip cycle if we already registered this center
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
    def cart_bf_num(self):
        return sum([shell.size() for shell in self.shells])

    def from_basis(self, name, shells_cls=None, **kwargs):
        from pysisyphus.wavefunction.Basis import shells_with_basis

        atoms, coords3d = self.atoms_coords3d
        if shells_cls is None:
            shells_cls = type(self)
        return shells_with_basis(
            atoms, coords3d, name=name, shells_cls=shells_cls, **kwargs
        )

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

    @staticmethod
    def from_pyscf_mol(mol):
        shells = list()
        for bas_id in range(mol.nbas):
            L = mol.bas_angular(bas_id)
            center = mol.bas_coord(bas_id)
            coeffs = mol.bas_ctr_coeff(bas_id).flatten()
            exps = mol.bas_exp(bas_id)
            assert (
                coeffs.size == exps.size
            ), "General contractions are not yet supported."
            center_ind = mol.bas_atom(bas_id)
            atom_symbol = mol.atom_symbol(center_ind)
            atomic_num = ATOMIC_NUMBERS[atom_symbol.lower()]
            shell = Shell(L, center, coeffs, exps, center_ind, atomic_num)
            shells.append(shell)
        return PySCFShells(shells)

    def center_shell_iter(self):
        sorted_shells = sorted(self.shells, key=lambda shell: shell.center_ind)
        return it.groupby(sorted_shells, key=lambda shell: shell.center_ind)

    def fragment_ao_map(self, fragments):
        frag_map = dict()
        for i, frag in enumerate(fragments):
            for center_ind in frag:
                frag_map[center_ind] = i

        frag_ao_map = dict()
        for i, ao_center in enumerate(self.sph_ao_centers):
            frag_ind = frag_map[ao_center]
            frag_ao_map.setdefault(frag_ind, list()).append(i)
        return frag_ao_map

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
        """Permutation matrix for spherical basis functions."""
        return sp.linalg.block_diag(*[self.sph_Ps[shell.L] for shell in self.shells])

    @property
    def P_cart(self):
        """Permutation matrix for Cartesian basis functions."""
        return sp.linalg.block_diag(*[self.cart_Ps[shell.L] for shell in self.shells])

    def eval(self, xyz, spherical=True):
        """Evaluate all basis functions at points xyz using generated code.

        A possibly more efficient approach is discussed in III C of
        https://doi.org/10.1063/1.469408.
        """
        if spherical:
            precontr = self.cart2sph_coeffs.T @ self.P_sph.T
        else:
            precontr = self.P_cart.T

        coords3d = self.coords3d
        all_vals = list()
        for center_ind, shells in self.center_shell_iter():
            center = coords3d[center_ind]
            Xa, Ya, Za = (xyz - center).T
            for shell in shells:
                L_tot, _, contr_coeffs, exponents = shell.as_tuple()
                # shape (nbfs in shell, number of primitives, number of points)
                vals = cart_gto(L_tot, exponents[:, None], Xa, Ya, Za)
                vals *= contr_coeffs[..., None]
                vals = vals.sum(axis=1)  # Contract primitives
                all_vals.append(vals)
        all_vals = np.concatenate(all_vals, axis=0).T
        all_vals = all_vals @ precontr
        return all_vals

    def eval_manual(self, xyz, spherical=False):
        """Evaluate all basis functions at points xyz using handwritten code."""
        all_vals = list()
        if spherical:
            precontr = self.cart2sph_coeffs.T @ self.P_sph.T
        else:
            precontr = self.P_cart.T

        for shell in self.shells:
            _, center, contr_coeffs, exponents = shell.as_tuple()
            cart_powers = shell.cart_powers
            vals = eval_shell(xyz, center, cart_powers, contr_coeffs, exponents)
            all_vals.append(vals)
        all_vals = np.concatenate(all_vals, axis=0).T
        all_vals = all_vals @ precontr
        return all_vals

    def get_1el_ints_cart(
        self,
        func,
        other=None,
        can_reorder=True,
        components=0,
        screen_func=None,
        **kwargs,
    ):
        shells_a = self
        symmetric = other is None
        shells_b = shells_a if symmetric else other
        cart_bf_num_a = shells_a.cart_bf_num
        cart_bf_num_b = shells_b.cart_bf_num

        # components = 0 indicates, that a plain 2d matrix is desired.
        if is_2d := (components == 0):
            components = 1

        # Preallocate empty matrices and directly assign the calculated values
        integrals = np.zeros((components, cart_bf_num_a, cart_bf_num_b))

        # Dummy screen function that always returns True
        if (not self.screen) or (screen_func is None):
            screen_func = lambda *args: True

        a_ind = 0
        b_skipped = 0
        for i, shell_a in enumerate(shells_a):
            La, A, dA, aa = shell_a.as_tuple()
            a_size = get_shell_shape(La)[0]
            a_slice = slice(a_ind, a_ind + a_size)
            b_ind = b_skipped
            if symmetric:
                b_skipped += a_size
            # When we don't deal with symmetric matrices we have to iterate over
            # all other basis functions in shells_b, so we reset i to 0.
            else:
                i = 0
            a_min_exp = aa.min()  # Minimum exponent used for screening
            for shell_b in shells_b[i:]:
                Lb, B, dB, bb = shell_b.as_tuple()
                b_min_exp = bb.min()  # Minimum exponent used for screening
                R_ab = np.linalg.norm(A - B)  # Shell center distance
                b_size = get_shell_shape(Lb)[0]
                if not screen_func(a_min_exp, b_min_exp, R_ab):
                    b_ind += b_size
                    continue
                shell_shape = (a_size, b_size)
                shape = (components, *shell_shape)
                ints_ = func(La, Lb, aa[:, None], A, bb[None, :], B, **kwargs)
                ints_ = ints_.reshape(*shape, len(aa), len(bb))  # 5d
                ints_ *= dA[None, :, None, :, None] * dB[None, None, :, None, :]
                ints_ = ints_.sum(axis=(3, 4))
                b_slice = slice(b_ind, b_ind + b_size)
                integrals[:, a_slice, b_slice] = ints_
                # Fill up the lower tringular part of the matrix
                if symmetric:
                    integrals[:, b_slice, a_slice] = np.transpose(ints_, axes=(0, 2, 1))
                b_ind += b_size
            a_ind += a_size

        # Return plain 2d array if components is set to 0, i.e., remove first axis.
        if is_2d:
            integrals = np.squeeze(integrals, axis=0)

        # Reordering will be disabled, when spherical integrals are desired. They
        # are reordered outside of this function. Reordering them already here
        # would mess up the results after the 2nd reordering.
        if can_reorder and self.ordering == "native":
            integrals = np.einsum(
                "ij,...jk,kl->...il",
                self.P_cart,
                integrals,
                shells_b.P_cart.T,
                optimize="greedy",
            )
        return integrals

    def get_1el_ints_sph(
        self, int_func=None, int_matrix=None, other=None, **kwargs
    ) -> NDArray:
        if int_matrix is None:
            # Disallow reordering of the Cartesian integrals
            int_matrix = self.get_1el_ints_cart(
                int_func, other=other, can_reorder=False, **kwargs
            )

        shells_b = self if other is None else other
        # Reordering (P) and Cartesian to spherical conversion (C).
        PC_a = self.reorder_c2s_coeffs
        PC_b = shells_b.reorder_c2s_coeffs
        int_matrix_sph = np.einsum(
            "ij,...jk,kl->...il", PC_a, int_matrix, PC_b.T, optimize="greedy"
        )
        return int_matrix_sph

    def get_2c2el_ints_cart(self):
        return self.get_1el_ints_cart(_2center2electron)

    def get_2c2el_ints_sph(self):
        return self.get_1el_ints_sph(_2center2electron)

    def get_3c2el_ints_cart(self, shells_aux, int_func=_3center2electron):
        shells_a = self
        cart_bf_num_a = shells_a.cart_bf_num
        cart_bf_num_aux = shells_aux.cart_bf_num

        integrals = np.zeros((cart_bf_num_a, cart_bf_num_a, cart_bf_num_aux))
        symmetric = True

        a_ind = 0
        b_skipped = 0
        for i, shell_a in enumerate(shells_a):
            La, A, dA, aa = shell_a.as_tuple()
            a_size = get_shell_shape(La)[0]
            a_slice = slice(a_ind, a_ind + a_size)
            b_ind = b_skipped
            if symmetric:
                b_skipped += a_size
            else:
                i = 0
            for shell_b in shells_a[i:]:
                Lb, B, dB, bb = shell_b.as_tuple()
                b_size = get_shell_shape(Lb)[0]
                b_slice = slice(b_ind, b_ind + b_size)
                dAB = (
                    dA[:, None, None, :, None, None] * dB[None, :, None, None, :, None]
                )
                c_ind = 0
                for shell_c in shells_aux:
                    Lc, C, dC, cc = shell_c.as_tuple()
                    c_size = get_shell_shape(Lc)[0]
                    c_slice = slice(c_ind, c_ind + c_size)
                    shape = (a_size, b_size, c_size)
                    ints_ = int_func(
                        La,
                        Lb,
                        Lc,
                        aa[:, None, None],
                        A,
                        bb[None, :, None],
                        B,
                        cc[None, None, :],
                        C,
                    )
                    ints_ = ints_.reshape(*shape, len(aa), len(bb), len(cc))  # 6D
                    ints_ *= dAB * dC[None, None, :, None, None, :]
                    ints_ = ints_.sum(axis=(3, 4, 5))
                    integrals[a_slice, b_slice, c_slice] = ints_
                    if symmetric:
                        integrals[b_slice, a_slice, c_slice] = np.transpose(
                            ints_, axes=(1, 0, 2)
                        )
                    c_ind += c_size
                b_ind += b_size
            a_ind += a_size
        return integrals

    def get_3c2el_ints_sph(self, shells_aux):
        int_matrix = self.get_3c2el_ints_cart(
            shells_aux, int_func=_3center2electron_sph
        )

        PC_a = self.reorder_c2s_coeffs
        PC_c = shells_aux.reorder_c2s_coeffs
        int_matrix_sph = np.einsum(
            "ij,jlm,nl,om->ino", PC_a, int_matrix, PC_a, PC_c, optimize="greedy"
        )
        return int_matrix_sph

    #####################
    # Overlap integrals #
    #####################

    @staticmethod
    def screen_S(a, b, R, k=10):
        """
        Returns False when distance R is below calculated threshold.
        Derived from the sympy code below.

        from sympy import *

        a, b, S, R = symbols("a b S R", real=True, positive=True)
        k = symbols("k", integer=True, positive=True)
        apb = a + b
        # 0 = S - 10**(-k)
        ovlp = (pi / apb)**Rational(3,2) * exp(-a*b / apb * R**2) - 10**(-k)
        print("ss-overlap:", ovlp)
        # Distance R, when ss-overlaps drops below 10**(-k)
        R = solve(ovlp, R)[0]
        print("R:", R)

        Parameters
        ----------
        a
            Minimum exponent in shell a.
        b
            Minimum exponent in shell b.
        R
            Real space distance of shells a and b.
        """
        return R < (
            sqrt(a + b)
            * sqrt(log(10**k * pi ** (3 / 2) / (a + b) ** (3 / 2)))
            / (sqrt(a) * sqrt(b))
        )

    def get_S_cart(self, other=None) -> NDArray:
        return self.get_1el_ints_cart(ovlp, other=other, screen_func=Shells.screen_S)

    def get_S_sph(self, other=None) -> NDArray:
        return self.get_1el_ints_sph(ovlp, other=other, screen_func=Shells.screen_S)

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

    def get_V_cart(self, coords3d=None, charges=None, **kwargs):
        # Coordinates and charges must be omitted/given together
        # Alternatively, this function could also take one array that
        # combines coords3c and charges.
        assert ((coords3d is None) and (charges is None)) or (
            (coords3d is not None) and (charges is not None)
        )
        if coords3d is None:
            atoms, coords3d = self.atoms_coords3d
            charges = nuc_charges_for_atoms(atoms)

        cart_bf_num = self.cart_bf_num
        V_nuc = np.zeros((cart_bf_num, cart_bf_num))
        # Loop over all centers and add their contributions
        for C, Z in zip(coords3d, charges):
            # -Z = -1 * Z, because electrons have negative charge.
            V_nuc += -Z * self.get_1el_ints_cart(
                coulomb,
                C=C,
                **kwargs,
            )
        return V_nuc

    def get_V_sph(self, coords3d=None, charges=None) -> NDArray:
        V_cart = self.get_V_cart(coords3d, charges, can_reorder=False)
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
        return self.get_1el_ints_cart(
            dipole, components=3, C=origin, screen_func=Shells.screen_S
        )

    def get_dipole_ints_sph(self, origin) -> NDArray:
        return self.get_1el_ints_sph(
            dipole, components=3, C=origin, screen_func=Shells.screen_S
        )

    ##################################################
    # Quadrupole moment integrals, diagonal elements #
    ##################################################

    def get_diag_quadrupole_ints_cart(self, origin):
        return self.get_1el_ints_cart(
            diag_quadrupole, components=3, C=origin, screen_func=Shells.screen_S
        )

    def get_diag_quadrupole_ints_sph(self, origin):
        return self.get_1el_ints_sph(
            diag_quadrupole, components=3, C=origin, screen_func=Shells.screen_S
        )

    ###############################
    # Quadrupole moment integrals #
    ###############################

    def get_quadrupole_ints_cart(self, origin):
        _ = self.get_1el_ints_cart(
            quadrupole, components=6, C=origin, screen_func=Shells.screen_S
        )
        shape = _.shape
        sym = np.zeros((3, 3, *shape[1:]))
        triu = np.triu_indices(3)
        triu1 = np.triu_indices(3, k=1)
        tril1 = np.tril_indices(3, k=-1)
        sym[triu] = _
        sym[tril1] = sym[triu1]
        return sym.reshape(3, 3, *shape[1:])

    def get_quadrupole_ints_sph(self, origin) -> NDArray:
        _ = self.get_1el_ints_sph(
            quadrupole, components=6, C=origin, screen_func=Shells.screen_S
        )
        shape = _.shape
        sym = np.zeros((3, 3, *shape[1:]))
        triu = np.triu_indices(3)
        triu1 = np.triu_indices(3, k=1)
        tril1 = np.tril_indices(3, k=-1)
        sym[triu] = _
        sym[tril1] = sym[triu1]
        return sym.reshape(3, 3, *shape[1:])

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
            [0, 0, 0, 0, 0, 0, -1],  # sign flip, as in the line above
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


class ORCAMoldenShells(Shells):
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
            [-1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, -1],
        ],
        4: [
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, -1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, -1, 0],
            [-1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -1],
        ],
    }


def pyscf_cart_order(l):
    order = list()
    for lx in reversed(range(l + 1)):
        for ly in reversed(range(l + 1 - lx)):
            lz = l - lx - ly
            order.append("x" * lx + "y" * ly + "z" * lz)
    return tuple(order)


class PySCFShells(Shells):

    """
    Cartesian bfs >= d angular momentum are not normalized!
    S_ref = mol.intor("int1e_ovlp_cart")
    N = 1 / np.diag(S_ref)**0.5
    ao *= N
    """

    cart_order = [pyscf_cart_order(l) for l in range(5)]

    sph_Ps = {
        0: [[1]],  # s
        1: [[1, 0, 0], [0, 0, 1], [0, 1, 0]],  # px py pz
        2: [
            [0, 0, 0, 0, 1],  # dxy
            [0, 0, 0, 1, 0],  # dyz
            [0, 0, 1, 0, 0],  # dz²
            [0, 1, 0, 0, 0],  # dxz
            [1, 0, 0, 0, 0],  # dx² - y²
        ],
        3: [
            [0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0],
        ],
        4: [
            [0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
    }
