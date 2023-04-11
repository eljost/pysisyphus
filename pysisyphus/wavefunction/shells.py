# [1] https://aip.scitation.org/doi/pdf/10.1063/5.0004046
#     Efficient implementation of the superposition of atomic potentials initial
#     guess for electronic structure calculations in Gaussian basis sets
#     Lehtola, Visscher, Engel,  2020

import itertools as it
from math import sqrt, log, pi
from pathlib import Path
import textwrap
from typing import List, Literal


from jinja2 import Template
from joblib import Memory
import numpy as np
from numpy.typing import NDArray
import scipy as sp


from pysisyphus.config import L_MAX, L_AUX_MAX
from pysisyphus.elem_data import (
    ATOMIC_NUMBERS,
    INV_ATOMIC_NUMBERS,
    nuc_charges_for_atoms,
)
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.linalg import multi_component_sym_mat
from pysisyphus.wavefunction.helpers import (
    canonical_order,
    get_l,
    get_shell_shape,
    L_MAP_INV,
    permut_for_order,
)

from pysisyphus.wavefunction.ints import (
    # boys,
    cart_gto3d,
    coulomb3d,
    diag_quadrupole3d,
    dipole3d,
    kinetic3d,
    ovlp3d,
    quadrupole3d,
    self_ovlp3d,
    _2center2el3d,
    _3center2el3d_sph,
)

from pysisyphus.wavefunction.cart2sph import cart2sph_coeffs


class Shell:
    def __init__(
        self,
        L: int,
        center: NDArray[float],
        coeffs: NDArray[float],
        exps: NDArray[float],
        center_ind: int,
        atomic_num=None,
        shell_index=None,
        index=None,
    ):
        self.L = get_l(L)
        self.center = np.array(center, dtype=float)  # (x, y, z), 1d array
        coeffs = np.array(coeffs, dtype=float)
        # Store original contraction coefficients
        self.coeffs_org = coeffs.copy()
        # Orbital exponents, 1d array
        self.exps = np.array(exps, dtype=float)
        assert self.coeffs_org.size == self.exps.size
        coeffs = np.array(coeffs)
        assert coeffs.size == self.exps.size
        self.center_ind = int(center_ind)
        if atomic_num is None:
            atomic_num = -1
        self.atomic_num = int(atomic_num)
        if shell_index is not None:
            self.shell_index = shell_index
        if index is not None:
            self.index = index

        # TODO: fix sympleints to NOT return 2d array
        self_overlaps = self_ovlp3d.self_ovlp3d[(self.L, self.L)](
            self.exps[:, None],
            coeffs[:, None],
            self.center,
            self.exps[None, :],
            coeffs[None, :],
            self.center,
        )[0]
        # Check, that all self overlaps in a shell are of the same value.
        np.testing.assert_allclose(
            self_overlaps - self_overlaps[0], np.zeros_like(self_overlaps), atol=1e-12
        )
        # Only use first self-overlap for normalization, as all values are equal.
        N = 1 / np.sqrt(self_overlaps[0])
        self.coeffs = N * coeffs

        try:
            self.atom = INV_ATOMIC_NUMBERS[self.atomic_num]
        except KeyError:
            self.atom = "X"
        self.atom = self.atom.lower()

    def as_tuple(self):
        return self.L, self.center, self.coeffs, self.exps, self.index, self.size

    def exps_coeffs_iter(self):
        return zip(self.exps, self.coeffs)

    @property
    def contr_depth(self):
        return self.coeffs.size

    @property
    def cart_powers(self):
        return np.array(canonical_order(self.L), dtype=int)

    @property
    def size(self):
        """Number of cartesian basis functions in the shell."""
        return self.cart_size

    @property
    def cart_size(self):
        """Number of cartesian basis functions in the shell."""
        L = self.L
        return (L + 1) * (L + 2) // 2

    @property
    def sph_size(self):
        """Number of cartesian basis functions in the shell."""
        return 2 * self.L + 1

    @property
    def index(self) -> int:
        return self._index

    @index.setter
    def index(self, index: int):
        self._index = int(index)

    @property
    def shell_index(self) -> int:
        return self._shell_index

    @shell_index.setter
    def shell_index(self, shell_index: int):
        self._shell_index = shell_index

    def to_pyscf_shell(self):
        key = L_MAP_INV[self.L]
        lines = [
            f"{self.atom.title()}    {key.upper()}",
        ]
        lines += [
            f"{exp_:> 18.10f}    {coeff:>18.10f}"
            for exp_, coeff in zip(self.exps, self.coeffs)
        ]
        return "\n".join(lines)

    def to_sap_shell(self):
        """Fix contraction coefficients, for use in SAP-initial guess calculation.

        See [1]."""
        # SAP shells must always be s-shells
        assert self.L == 0
        self.coeffs = self.coeffs_org
        coeffs_sum = self.coeffs.sum()
        # Check sum rule; sum should equal -Z.
        assert (coeffs_sum + ATOMIC_NUMBERS[self.atom.lower()]) <= 1e-5
        # Fix normalization. The primitive Gaussians are normalized inside
        # the integral functions, so we pre-divide the contraction coefficients
        # by the normalization factor.
        N = (2 * self.exps / np.pi) ** 0.75
        self.coeffs /= N
        # The potential fits were carried out for functions of the form:
        #   g_p(r) = (α_p / π)**1.5 * exp(-α_p / r) .
        # So we multiply the prefactor onto the contraction coefficients.
        # See [1], p. 152, just above Eq. (2).
        self.coeffs *= (self.exps / np.pi) ** 1.5

    def __str__(self):
        try:
            center_str = f", at atom {self.center_ind}"
        except AttributeError:
            center_str = ""
        return f"Shell(L={self.L}, {self.contr_depth} pGTO{center_str})"

    def __repr__(self):
        return self.__str__()


Ordering = Literal["native", "pysis"]


class Shells:
    sph_Ps = {l: np.eye(2 * l + 1) for l in range(max(L_MAX, L_AUX_MAX) + 1)}

    def __init__(
        self,
        shells: List[Shell],
        screen: bool = False,
        cache: bool = False,
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

        # Now that we have all Shell objects, we can set their starting indices
        index = 0
        shell_index = 0
        for shell in self.shells:
            shell.shell_index = shell_index
            shell.index = index
            index += shell.size
            shell_index += 1

        # Try to construct Cartesian permutation matrix from cart_order, if defnied.
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

    def as_tuple(self):
        # Ls, centers, contr_coeffs, exponents, indices, sizes
        return zip(*[shell.as_tuple() for shell in self.shells])

    @property
    def cart_size(self):
        """Number of cartesian basis functions."""
        return sum([shell.cart_size for shell in self.shells])

    @property
    def sph_size(self):
        """Number of spherical basis."""
        return sum([shell.sph_size for shell in self.shells])

    def as_arrays(self, fortran=False):
        """Similar to the approach in libcint."""
        # bas_centers
        # One entry per shell, integer array.
        #   center_ind, atomic number, pointer to center coordinates in bas_data (3 integers)
        bas_centers = list()
        # bas_spec
        # One entry per shell, integer array.
        #   shell_ind, total angmom, N_pgto, N_cgto, \
        #   pointer to contraction coefficients and exponents in bas_data \
        #   (2*N_pgto floats)
        bas_spec = list()
        # bas_data
        # Float array. Starts with 3 * N_centers floats, containing the center
        # coordinates. Continues with alternating contraction coefficient and
        # orbital exponent data.
        bas_data = list()
        pointer = 1 if fortran else 0

        shells = self.shells
        # Store center data, i.e., where the basis functions are located.
        for shell in shells:
            center_ind = shell.center_ind
            bas_data.extend(shell.center)  # Store coordinates
            atomic_num = shell.atomic_num
            if atomic_num is None:
                atomic_num = -1
            bas_centers.append((center_ind, atomic_num, pointer))
            pointer += 3

        # Store contraction coefficients and primitive exponents.
        for shell_ind, shell in enumerate(self.shells):
            # L, center, coeffs, exponents, index, size = shell.as_tuple()
            L, _, _, exponents, _, _ = shell.as_tuple()
            contr_coeffs = shell.coeffs_org
            nprim = len(contr_coeffs)
            # ncontr is hardcoded to 1 for now, as there are no special functions
            # to handle general contractions.
            ncontr = 1
            bas_spec.append((shell_ind, L, nprim, ncontr, pointer))
            bas_data.extend(contr_coeffs)
            pointer += nprim
            bas_data.extend(exponents)
            pointer += nprim

        bas_centers = np.array(bas_centers, dtype=int)
        bas_spec = np.array(bas_spec, dtype=int)
        bas_data = np.array(bas_data, dtype=float)
        return bas_centers, bas_spec, bas_data

    def as_sympleints_arrays(self):
        shells = self.shells
        # center_ind, atomic_num, L, nprims
        centers = np.zeros((len(shells), 3))
        shell_data = list()
        coefficients = list()
        exponents = list()
        for i, shell in enumerate(shells):
            coeffs = shell.coeffs_org
            coefficients.extend(coeffs.tolist())
            exponents.extend(shell.exps.tolist())
            nprims = len(coeffs)

            centers[i] = shell.center
            atomic_num = shell.atomic_num
            if atomic_num is None:
                atomic_num = -1
            shell_data.append((shell.center_ind, atomic_num, shell.L, nprims))
        shell_data = np.array(shell_data, dtype=np.int32)
        coefficients = np.array(coefficients)
        exponents = np.array(exponents)
        return shell_data, centers, coefficients, exponents

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
        return sum([shell.size for shell in self.shells])

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
        from pysisyphus.io.aomix import shells_from_aomix

        shells = shells_from_aomix(text)
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

    def to_pyscf_mol(self):
        # TODO: this would be better suited as a method of Wavefunction,
        # as pyscf Moles must have sensible spin & charge etc.
        # TODO: currently, this does not support different basis sets
        # for the same element.
        try:
            from pyscf import gto
        except ModuleNotFoundError:
            return None

        unique_atoms = set()
        basis = dict()
        atoms_coords = list()
        for _, center_shells in self.center_shell_iter():
            center_shells = list(center_shells)
            shell0 = center_shells[0]
            atom = shell0.atom
            x, y, z = shell0.center
            atoms_coords.append(f"{atom} {x} {y} {z}")
            if atom not in unique_atoms:
                unique_atoms.add(atom)
            else:
                continue
            atom_shells = "\n".join([shell.to_pyscf_shell() for shell in center_shells])
            basis[atom] = gto.basis.parse(atom_shells)

        mol = gto.Mole()
        mol.atom = "; ".join(atoms_coords)
        mol.unit = "Bohr"
        mol.basis = basis
        mol.build()
        return mol

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
        try:
            P_cart = sp.linalg.block_diag(
                *[self.cart_Ps[shell.L] for shell in self.shells]
            )
        except AttributeError:
            P_cart = np.eye(self.cart_size)
        return P_cart

    def eval(self, xyz, spherical=True):
        """Evaluate all basis functions at points xyz using generated code.

        A possibly more efficient approach is discussed in III C of
        https://doi.org/10.1063/1.469408.
        """
        if spherical:
            precontr = self.cart2sph_coeffs.T @ self.P_sph.T
        else:
            precontr = self.P_cart.T
        ncbfs = precontr.shape[0]
        npoints = len(xyz)

        coords3d = self.coords3d
        grid_vals = np.zeros((npoints, ncbfs))
        for center_ind, shells in self.center_shell_iter():
            center = coords3d[center_ind]
            # Create list from the shells iterator, otherwise it will be empty
            # after the first pass over all points.
            shells = list(shells)
            for i, (X, Y, Z) in enumerate(xyz - center):
                for shell in shells:
                    La, A, da, ax, a_ind, a_size = shell.as_tuple()
                    grid_vals[i, a_ind : a_ind + a_size] = cart_gto3d.cart_gto3d[(La,)](
                        ax, da, X, Y, Z
                    )
        return grid_vals @ precontr

    def get_1el_ints_cart(
        self,
        func_dict,
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

        for i, shell_a in enumerate(shells_a):
            La, A, da, ax, a_ind, a_size = shell_a.as_tuple()
            a_slice = slice(a_ind, a_ind + a_size)
            # When we don't deal with symmetric matrices we have to iterate over
            # all other basis functions in shells_b, so we reset i to 0.
            if not symmetric:
                i = 0
            a_min_exp = ax.min()  # Minimum exponent used for screening
            for shell_b in shells_b[i:]:
                Lb, B, db, bx, b_ind, b_size = shell_b.as_tuple()
                b_slice = slice(b_ind, b_ind + b_size)
                # Screen shell pair
                b_min_exp = bx.min()  # Minimum exponent used for screening
                R_ab = np.linalg.norm(A - B)  # Shell center distance
                b_size = get_shell_shape(Lb)[0]
                if not screen_func(a_min_exp, b_min_exp, R_ab):
                    continue
                integrals[:, a_slice, b_slice] = func_dict[(La, Lb)](
                    ax[:, None],
                    da[:, None],
                    A,
                    bx[None, :],
                    db[None, :],
                    B,
                    **kwargs,
                )
                # Fill up the lower tringular part of the matrix in the symmetric
                # case.
                # TODO: this may be (?) better done outside of this loop, after
                # everything is calculated...
                if symmetric:
                    integrals[:, b_slice, a_slice] = np.transpose(
                        integrals[:, a_slice, b_slice], axes=(0, 2, 1)
                    )

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
        return self.get_1el_ints_cart(_2center2el3d._2center2el3d)

    def get_2c2el_ints_sph(self):
        return self.get_1el_ints_sph(_2center2el3d._2center2el3d)

    def get_3c2el_ints_cart(self, shells_aux):
        """Cartesian 3-center-2-electron integrals.

        DO NOT USE THESE INTEGRALS AS THEY ARE RETURNED FROM THIS METHOD.
        These integrals utilize recurrence relations that are only valid,
        when the resulting Cartesian integrals are transformed into spherical
        integrals.

        Contrary to the general function 'get_1el_ints_cart', that supports
        different 'func_dict' arguments and cross-integrals between two
        different shells this function is less general. This function is
        restricted to '_3center2el_sph' and always uses 'self.shells' as well as
        'self.shells_aux'.
        """
        shells_a = self
        cart_bf_num_a = shells_a.cart_bf_num
        cart_bf_num_aux = shells_aux.cart_bf_num

        integrals = np.zeros((cart_bf_num_a, cart_bf_num_a, cart_bf_num_aux))
        func_dict = _3center2el3d_sph._3center2el3d_sph

        for i, shell_a in enumerate(shells_a):
            La, A, da, ax, a_ind, a_size = shell_a.as_tuple()
            a_slice = slice(a_ind, a_ind + a_size)
            # As noted in the docstring, we iterate over pairs of self.shells (shells_a)
            for shell_b in shells_a[i:]:
                Lb, B, db, bx, b_ind, b_size = shell_b.as_tuple()
                b_slice = slice(b_ind, b_ind + b_size)
                for shell_c in shells_aux:
                    Lc, C, dc, cx, c_ind, c_size = shell_c.as_tuple()
                    c_slice = slice(c_ind, c_ind + c_size)
                    integrals[a_slice, b_slice, c_slice] = func_dict[(La, Lb, Lc)](
                        ax[:, None, None],
                        da[:, None, None],
                        A,
                        bx[None, :, None],
                        db[None, :, None],
                        B,
                        cx[None, None, :],
                        dc[None, None, :],
                        C,
                    )
                    integrals[b_slice, a_slice, c_slice] = np.transpose(
                        integrals[a_slice, b_slice, c_slice], axes=(1, 0, 2)
                    )
        return integrals

    def get_3c2el_ints_sph(self, shells_aux):
        int_matrix = self.get_3c2el_ints_cart(shells_aux)

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
        # return self.get_1el_ints_cart(ovlp, other=other, screen_func=Shells.screen_S)
        return self.get_1el_ints_cart(
            ovlp3d.ovlp3d, other=other, screen_func=Shells.screen_S
        )

    def get_S_sph(self, other=None) -> NDArray:
        return self.get_1el_ints_sph(
            ovlp3d.ovlp3d, other=other, screen_func=Shells.screen_S
        )

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
        # return self.get_1el_ints_cart(kinetic, other=other)
        return self.get_1el_ints_cart(kinetic3d.kinetic3d, other=other)

    def get_T_sph(self, other=None) -> NDArray:
        return self.get_1el_ints_sph(kinetic3d.kinetic3d, other=other)

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
        """Nuclear attraction integrals.

        Coordinates and charges must be either be omitted or given together.
        Alternatively, this function could also take one array that combines
        coords3c and charges (not yet implemented).
        """
        assert ((coords3d is None) and (charges is None)) or (
            (coords3d is not None) and (charges is not None)
        )
        if coords3d is None:
            atoms, coords3d = self.atoms_coords3d
            charges = nuc_charges_for_atoms(atoms)

        # def boys_1c1el(n_max, ax, A, bx, B, C):
        # ax = ax[:, None]
        # bx = bx[None, :]
        # p = ax + bx
        # P = -(ax * A + bx * B) / p
        # R_PC2 = ((P - C)**2).sum()
        # pR_PC2 = p * R_PC2
        # return np.array([boys.boys(n, pR_PC2) for n in range(n_max+1)])

        cart_bf_num = self.cart_bf_num
        V_nuc = np.zeros((cart_bf_num, cart_bf_num))
        # Loop over all centers and add their contributions
        for C, Z in zip(coords3d, charges):
            # -Z = -1 * Z, because electrons have negative charge.
            V_nuc += -Z * self.get_1el_ints_cart(
                coulomb3d.coulomb3d,
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
            dipole3d.dipole3d, components=3, R=origin, screen_func=Shells.screen_S
        )

    def get_dipole_ints_sph(self, origin) -> NDArray:
        return self.get_1el_ints_sph(
            dipole3d.dipole3d, components=3, R=origin, screen_func=Shells.screen_S
        )

    ##################################################
    # Quadrupole moment integrals, diagonal elements #
    ##################################################

    def get_diag_quadrupole_ints_cart(self, origin):
        return self.get_1el_ints_cart(
            diag_quadrupole3d.diag_quadrupole3d,
            components=3,
            R=origin,
            screen_func=Shells.screen_S,
        )

    def get_diag_quadrupole_ints_sph(self, origin):
        return self.get_1el_ints_sph(
            diag_quadrupole3d.diag_quadrupole3d,
            components=3,
            R=origin,
            screen_func=Shells.screen_S,
        )

    ###############################
    # Quadrupole moment integrals #
    ###############################

    def get_quadrupole_ints_cart(self, origin):
        ints_flat = self.get_1el_ints_cart(
            quadrupole3d.quadrupole3d,
            components=6,
            R=origin,
            screen_func=Shells.screen_S,
        )
        return multi_component_sym_mat(ints_flat, 3)

    def get_quadrupole_ints_sph(self, origin) -> NDArray:
        ints_flat = self.get_1el_ints_sph(
            quadrupole3d.quadrupole3d,
            components=6,
            R=origin,
            screen_func=Shells.screen_S,
        )
        return multi_component_sym_mat(ints_flat, 3)

    def to_sap_shells(self):
        for shell in self.shells:
            shell.to_sap_shell()

    def __str__(self):
        return f"{self.__class__.__name__}({len(self.shells)} shells, ordering={self.ordering})"


class AOMixShells(Shells):
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
