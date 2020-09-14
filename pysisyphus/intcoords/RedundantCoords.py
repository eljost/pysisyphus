# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844

import itertools as it
import logging

import numpy as np

from pysisyphus.helpers_pure import remove_duplicates
from pysisyphus.intcoords import Bend, LinearBend, Stretch, Torsion
from pysisyphus.intcoords.backconversion import transform_int_step
from pysisyphus.intcoords.derivatives import d2q_b, d2q_a, d2q_d
from pysisyphus.intcoords.eval import (
    eval_primitives,
    check_primitives,
    augment_primitives,
    PrimInternal,
)
from pysisyphus.intcoords.setup import setup_redundant, get_primitives, valid_bend, valid_dihedral


class RedundantCoords:
    def __init__(
        self,
        atoms,
        cart_coords,
        bond_factor=1.3,
        prim_indices=None,
        define_prims=None,
        bonds_only=False,
        check_bends=True,
        check_dihedrals=False,
        bend_min_deg=15,
        bend_max_deg=180,
        lb_min_deg=None,
        make_complement=True,
    ):
        self.atoms = atoms
        self.cart_coords = cart_coords
        self.bond_factor = bond_factor
        self.define_prims = define_prims
        self.bonds_only = bonds_only
        self.check_bends = check_bends
        self.check_dihedrals = check_dihedrals
        self.bend_min_deg = bend_min_deg
        self.bend_max_deg = bend_max_deg
        self.lb_min_deg = lb_min_deg
        self.make_complement = make_complement

        self._B_prim = None
        # Lists for the other types of primitives will be created afterwards.
        # Linear bends may have been disabled, so we create the list here.
        self.linear_bend_indices = list()
        self.logger = logging.getLogger("internal_coords")

        # Set up primitive indices
        if prim_indices is None:
            self.set_primitive_indices(
                self.atoms,
                self.coords3d,
                min_deg=self.bend_min_deg,
                max_deg=self.bend_max_deg,
                define_prims=self.define_prims,
            )
        else:
            to_arr = lambda _: np.array(list(_), dtype=int)
            bonds, bends, dihedrals = prim_indices
            # We accept all bond indices. What could possibly go wrong?! :)
            self.bond_indices = to_arr(bonds)
            valid_bends = [inds for inds in bends if valid_bend(self.coords3d, inds,
                self.bend_min_deg, self.bend_max_deg)]
            self.bending_indices = to_arr(valid_bends)
            valid_dihedrals = [
                inds for inds in dihedrals if valid_dihedral(self.coords3d, inds)
            ]
            self.dihedral_indices = to_arr(valid_dihedrals)

        if self.bonds_only:
            self.bending_indices = list()
            self.dihedral_indices = list()

        self.primitives = get_primitives(
            self.coords3d,
            self.bond_indices,
            self.bending_indices,
            self.linear_bend_indices,
            self.dihedral_indices,
            make_complement=self.make_complement,
            logger=self.logger,
        )

        self._prim_internals = self.eval(self.coords3d)
        self._prim_coords = np.array([prim_int.val for prim_int in self._prim_internals])

        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        self._bonds_slice = slice(bonds)
        self._bends_slice = slice(bonds, bonds+bends)
        self._dihedrals_slice = slice(bonds+bends, bonds+bends+dihedrals)

    def log(self, message):
        self.logger.debug(message)

    @property
    def cart_coords(self):
        return self._cart_coords

    @cart_coords.setter
    def cart_coords(self, cart_coords):
        self._cart_coords = cart_coords
        self._B_prim = None
        self._prim_coords = None

    @property
    def primitives(self):
        return self._primitives

    @primitives.setter
    def primitives(self, primitives):
        self._primitives = primitives

    @property
    def coords3d(self):
        return self.cart_coords.reshape(-1, 3)

    @property
    def prim_indices(self):
        return [self.bond_indices, self.bending_indices, self.dihedral_indices]

    @property
    def prim_indices_set(self):
        return set([tuple(prim_ind) for prim_ind in it.chain(*self.prim_indices)])

    @property
    def prim_internals(self):
        if self._prim_internals is None:
            self._prim_internals = self.eval(self.coords3d)
        return self._prim_internals

    @prim_internals.setter
    def prim_internals(self, prim_internals):
        self._prim_internals = prim_internals

    @property
    def prim_coords(self):
        return np.array([prim_int.val for prim_int in self.prim_internals])

    def return_inds(self, slice_):
        return np.array(
            [prim_int.indices for prim_int in self.prim_internals[slice_]]
        )

    @property
    def bonds(self):
        return self.prim_internals[self._bonds_slice]

    @property
    def bends(self):
        return self.prim_internals[self._bends_slice]

    @property
    def dihedrals(self):
        return self.prim_internals[self._dihedrals_slice]

    @property
    def coords(self):
        return self.prim_coords

    # @property
    # def coord_indices(self):
        # ic_ind_tuples = [tuple(prim.indices) for prim in self._primitives]
        # return {ic_inds: i for i, ic_inds in enumerate(ic_ind_tuples)}

    @property
    def dihed_start(self):
        return len(self.bond_indices) + len(self.bending_indices)

    def get_index_of_prim_coord(self, prim_ind):
        """Index of primitive internal for the given atom indices."""
        prim_ind_set = set(prim_ind)
        for i, prim in enumerate(self.primitives):
            if set(prim.indices) == prim_ind_set:
                return i
        self.log(f"Primitive internal with indices {prim_ind} " "is not defined!")
        return None

    @property
    def B_prim(self):
        """Wilson B-Matrix"""
        if self._B_prim is None:
            self._B_prim = np.array([prim_int.grad for prim_int in self.prim_internals])

        return self._B_prim

    @property
    def B(self):
        """Wilson B-Matrix"""
        return self.B_prim

    @property
    def Bt_inv_prim(self):
        """Transposed generalized inverse of the primitive Wilson B-Matrix."""
        B = self.B_prim
        return np.linalg.pinv(B.dot(B.T)).dot(B)

    @property
    def Bt_inv(self):
        """Transposed generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return np.linalg.pinv(B.dot(B.T)).dot(B)

    @property
    def B_inv_prim(self):
        """Generalized inverse of the primitive Wilson B-Matrix."""
        B = self.B_prim
        return B.T.dot(np.linalg.pinv(B.dot(B.T)))


    @property
    def B_inv(self):
        """Generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return B.T.dot(np.linalg.pinv(B.dot(B.T)))

    @property
    def P(self):
        """Projection matrix onto B. See [1] Eq. (4)."""
        return self.B.dot(self.B_inv)

    def transform_forces(self, cart_forces):
        """Combination of Eq. (9) and (11) in [1]."""
        return self.Bt_inv.dot(cart_forces)

    def get_K_matrix(self, int_gradient=None):
        if int_gradient is not None:
            assert len(int_gradient) == len(self._primitives)
        size_ = self.cart_coords.size
        if int_gradient is None:
            return np.zeros((size_, size_))

        dg_funcs = {
            2: d2q_b,
            # Todo: handle linear bend
            3: d2q_a,
            4: d2q_d,
        }

        def grad_deriv_wrapper(inds):
            coords_flat = self.coords3d[inds].flatten()
            dgrad = dg_funcs[len(inds)](*coords_flat)
            return dgrad

        K_flat = np.zeros(size_ * size_)
        for primitive, int_grad_item in zip(self.primitives, int_gradient):
            # Contract with gradient
            try:
                dg = int_grad_item * grad_deriv_wrapper(primitive.indices)
            except (ValueError, ZeroDivisionError) as err:
                self.log(
                    "Error in calculation of 2nd derivative of primitive "
                    f"internal {primitive.indices}."
                )
                continue
            # Depending on the type of internal coordinate dg is a flat array
            # of size 36 (stretch), 81 (bend) or 144 (torsion).
            #
            # An internal coordinate contributes to an element K[j, k] of the
            # K matrix if the cartesian coordinate indices j and k belong to an
            # atom that contributes to the respective internal coordinate.
            #
            # As for now we build up the K matrix as flat array. To add the dg
            # entries at the appropriate places in K_flat we have to calculate
            # the corresponding flat indices of dg in K_flat.
            cart_inds = list(it.chain(*[range(3 * i, 3 * i + 3) for i in primitive.indices]))
            flat_inds = [
                row * size_ + col for row, col in it.product(cart_inds, cart_inds)
            ]
            K_flat[flat_inds] += dg
        K = K_flat.reshape(size_, size_)
        return K

    def transform_hessian(self, cart_hessian, int_gradient=None):
        """Transform Cartesian Hessian to internal coordinates."""
        if int_gradient is None:
            self.log(
                "Supplied 'int_gradient' is None. K matrix will be zero, "
                "so derivatives of the Wilson-B-matrix are neglected in "
                "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv_prim.dot(cart_hessian - K).dot(self.B_inv_prim)

    def backtransform_hessian(self, redund_hessian, int_gradient=None):
        """Transform Hessian in internal coordinates to Cartesians."""
        if int_gradient is None:
            self.log(
                "Supplied 'int_gradient' is None. K matrix will be zero, "
                "so derivatives of the Wilson-B-matrix are neglected in "
                "the hessian transformation."
            )
        K = self.get_K_matrix(int_gradient)
        return self.B.T.dot(redund_hessian).dot(self.B) + K

    def project_hessian(self, H, shift=1000):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return P.dot(H).dot(P) + shift * (np.eye(P.shape[0]) - P)

    def project_vector(self, vector):
        """Project supplied vector onto range of B."""
        return self.P.dot(vector)

    def set_primitive_indices(
        self, atoms, coords3d, min_deg, max_deg, define_prims=None
    ):
        coord_info = setup_redundant(
            atoms,
            coords3d,
            factor=self.bond_factor,
            define_prims=define_prims,
            min_deg=min_deg,
            max_deg=max_deg,
            lb_min_deg=self.lb_min_deg,
            logger=self.logger,
        )

        all_bonds = (
            coord_info.bonds + coord_info.hydrogen_bonds + coord_info.interfrag_bonds
        )
        all_bonds = remove_duplicates(all_bonds)

        # Set primitive indices
        self.bond_indices = all_bonds
        self.bending_indices = coord_info.bends
        self.linear_bend_indices = coord_info.linear_bends
        self.dihedral_indices = coord_info.dihedrals

        self.hydrogen_bond_indices = coord_info.hydrogen_bonds
        self.fragments = coord_info.fragments
        self.cbm = coord_info.cbm
        self.cdm = coord_info.cdm

        # TODO: primitives are not yet defined
        # missing_prims, kappa = check_primitives(
            # self.coords3d,
            # self.primitives,
            # logger=self.logger,
        # )

    # def set_prim_internals(self, prim_internals):
        # self._prim_internals = prim_internals
        # self._prim_coords = np.array([prim.val for prim in self._prim_internals])

    def eval(self, coords3d, attr=None):
        prim_internals = eval_primitives(coords3d, self.primitives)

        if attr is not None:
            return np.array([
                getattr(prim_internal, attr) for prim_internal in prim_internals
            ])

        return prim_internals

    def transform_int_step(self, int_step, pure=False):
        new_prim_internals, cart_step, failed = transform_int_step(
            int_step,
            self.cart_coords,
            self.prim_coords,
            self.B_prim,
            self.primitives,
            logger=self.logger,
        )
        if not pure:
            self.prim_internals = new_prim_internals
        return cart_step

    def __str__(self):
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"
