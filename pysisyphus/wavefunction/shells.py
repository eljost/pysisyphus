from jinja2 import Template
import itertools as it
import textwrap
from typing import Tuple


import numpy as np
from numpy.typing import NDArray
import scipy as sp
from scipy.special import factorial2


from pysisyphus.config import L_MAX
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

from pysisyphus.wavefunction import (
    coulomb3d,
    diag_quadrupole3d,
    dipole3d,
    gto3d,
    kinetic3d,
    ovlp3d,
    quadrupole3d,
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
        self.coeffs = np.atleast_2d(coeffs)
        self.exps = np.array(exps, dtype=float)  # Orbital exponents, 1d array
        assert self.coeffs.shape[1] == self.exps.size
        self.center_ind = center_ind
        self.atomic_num = atomic_num

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


Ls = list(range(L_MAX + 1))


def get_map(module, func_base_name, Ls_num=2):
    """Return dict that holds the different integrals functions."""
    func_map = dict()
    for ls in it.product(*[Ls for _ in range(Ls_num)]):
        ls_str = "".join([str(l) for l in ls])
        func_map[ls] = getattr(module, f"{func_base_name}_{ls_str}")
    return func_map


CGTOmap = get_map(gto3d, "cart_gto3d", Ls_num=1)  # Cartesian GTO shells
Smap = get_map(ovlp3d, "ovlp3d")  # Overlap integrals
Tmap = get_map(kinetic3d, "kinetic3d")  # Kinetic energy integrals
Vmap = get_map(coulomb3d, "coulomb3d")  # 1el Coulomb integrals
DPMmap = get_map(dipole3d, "dipole3d")  # Dipole moments integrals
QPMmap = get_map(quadrupole3d, "quadrupole3d")  # Quadrupole moments integrals
DQPMmap = get_map(
    diag_quadrupole3d, "diag_quadrupole3d"
)  # Diagonal quadrupole moments integrals


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

        self.atoms, self.coords3d = self.atoms_coords3d

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

    @staticmethod
    def from_pyscf_mol(mol):
        shells = list()
        for bas_id in range(mol.nbas):
            L = mol.bas_angular(bas_id)
            center = mol.bas_coord(bas_id)
            coeffs = mol.bas_ctr_coeff(bas_id).flatten()
            exps = mol.bas_exp(bas_id)
            center_ind = mol.bas_atom(bas_id)
            atom_symbol = mol.atom_symbol(center_ind)
            atomic_num = ATOMIC_NUMBERS[atom_symbol.lower()]
            shell = Shell(L, center, coeffs, exps, center_ind, atomic_num)
            shells.append(shell)
        return PySCFShells(shells)

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
        self, int_func, other=None, add_args=None, can_reorder: bool = True
    ) -> NDArray:
        shells_a = self
        shells_b = shells_a if other is None else other

        if add_args is None:
            add_args = {}

        rows = list()
        for shell_a in shells_a.shells:
            La, A, dA, aa = shell_a.as_tuple()
            row = list()
            for shell_b in shells_b.shells:
                Lb, B, dB, bb = shell_b.as_tuple()
                shape = get_shell_shape(La, Lb)
                # Integral over primitives, shape: (product(*shape), prims_a, prims_b)
                pints = int_func(La, Lb, aa[:, None], A, bb[None, :], B, **add_args)
                pints = pints.reshape(*shape, len(aa), len(bb))
                pints *= dA[:, None, :, None] * dB[None, :, None, :]
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
        cart_bf_num = self.cart_bf_num
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

    def get_multipole_ints_cart(self, shells_a, shells_b, func, components, **kwargs):
        cart_bf_num = self.cart_bf_num
        # Preallocate empty matrices and directly assign the calculated values
        integrals = np.zeros((components, cart_bf_num, cart_bf_num))
        a_ind = 0
        for shell_a in shells_a.shells:
            La, A, dA, aa = shell_a.as_tuple()
            b_ind = 0
            for shell_b in shells_b.shells:
                Lb, B, dB, bb = shell_b.as_tuple()
                shell_shape = get_shell_shape(La, Lb)
                shape = (components, *shell_shape)
                a_size, b_size = shell_shape
                qp = func(La, Lb, aa[:, None], A, bb[None, :], B, **kwargs)
                qp = qp.reshape(*shape, len(aa), len(bb))
                qp *= dA[None, :, None, :, None] * dB[None, None, :, None, :]
                qp = qp.sum(axis=(3, 4)).reshape(shape)
                integrals[:, a_ind : a_ind + a_size, b_ind : b_ind + b_size] = qp
                b_ind += (Lb + 1) * (Lb + 2) // 2
            a_ind += (La + 1) * (La + 2) // 2
        return integrals

    def get_multipole_ints_sph(self, shells_a, shells_b, func, **kwargs) -> NDArray:
        cart_ints = self.get_multipole_ints_cart(shells_a, shells_b, func, **kwargs)
        C_a = shells_a.cart2sph_coeffs
        C_b = shells_b.cart2sph_coeffs
        sph_ints = list()
        for ci in cart_ints:
            ti = C_a.dot(ci).dot(C_b.T)  # Cartsian -> Spherical conversion
            if self.ordering == "native":
                ti = shells_a.P_sph.dot(ti).dot(shells_b.P_sph.T)  # Reorder
            sph_ints.append(ti)
        return np.array(sph_ints)

    ###########################
    # Dipole moment integrals #
    ###########################

    def get_dipole_ints_cart(self, origin):
        return self.get_multipole_ints_cart(self, self, dipole, components=3, C=origin)

    def get_dipole_ints_sph(self, origin) -> NDArray:
        return self.get_multipole_ints_sph(self, self, dipole, components=3, C=origin)

    ##################################################
    # Quadrupole moment integrals, diagonal elements #
    ##################################################

    def get_diag_quadrupole_ints_cart(self, origin):
        return self.get_multipole_ints_cart(
            self, self, diag_quadrupole, components=3, C=origin
        )

    def get_diag_quadrupole_ints_sph(self, origin):
        return self.get_multipole_ints_sph(
            self, self, diag_quadrupole, components=3, C=origin
        )

    ###############################
    # Quadrupole moment integrals #
    ###############################

    def get_quadrupole_ints_cart(self, origin):
        _ = self.get_multipole_ints_cart(self, self, quadrupole, components=6, C=origin)
        shape = _.shape
        sym = np.zeros((3, 3, *shape[1:]))
        triu = np.triu_indices(3)
        triu1 = np.triu_indices(3, k=1)
        tril1 = np.tril_indices(3, k=-1)
        sym[triu] = _
        sym[tril1] = sym[triu1]
        return sym.reshape(3, 3, *shape[1:])

    def get_quadrupole_ints_sph(self, origin) -> NDArray:
        _ = self.get_multipole_ints_sph(self, self, quadrupole, components=6, C=origin)
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
