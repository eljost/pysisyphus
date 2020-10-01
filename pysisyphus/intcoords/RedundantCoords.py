# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844

import math
import itertools as it
import logging

import numpy as np

from pysisyphus.intcoords import Stretch, Torsion
from pysisyphus.intcoords.update import transform_int_step
from pysisyphus.intcoords.eval import (
    eval_primitives,
    check_primitives,
)
from pysisyphus.intcoords.setup import setup_redundant, get_primitives, PrimTypes
from pysisyphus.intcoords.valid import check_typed_prims


class RedundantCoords:
    def __init__(
        self,
        atoms,
        coords3d,
        bond_factor=1.3,
        typed_prims=None,
        define_prims=None,
        bonds_only=False,
        check_bends=True,
        rebuild=True,
        bend_min_deg=15,
        dihed_max_deg=175.0,
        lb_min_deg=175.0,
        weighted=False,
        min_weight=0.3,
        rcond=1e-8,
    ):
        self.atoms = atoms
        self.coords3d = np.reshape(coords3d, (-1, 3)).copy()
        self.bond_factor = bond_factor
        self.define_prims = define_prims
        self.bonds_only = bonds_only
        self.check_bends = check_bends
        self.rebuild = rebuild
        self.bend_min_deg = bend_min_deg
        self.dihed_max_deg = dihed_max_deg
        self.lb_min_deg = lb_min_deg
        self.weighted = weighted
        self.min_weight = float(min_weight)
        assert self.min_weight > 0.0, "min_weight must be a positive rational!"
        self.rcond = rcond

        self._B_prim = None
        # Lists for the other types of primitives will be created afterwards.
        # Linear bends may have been disabled, so we create the list here.
        self.linear_bend_indices = list()
        self.logger = logging.getLogger("internal_coords")

        if self.weighted:
            self.log(
                "Coordinate weighting requested, min_weight="
                f"{self.min_weight:.2f}. Calculating bond factor."
            )
            # Screening function is
            #   ρ(d) = exp(-(d/sum_cov_rad - 1)
            #
            # Swart proposed a min_weight of ρ(d) = 0.3. With this we can
            # calculate the appropriate factor for the bond detection.
            # d = (1 - ln(0.3)) * sum_cov_rad
            # bond_factor = (1 - ln(0.3)) ≈ 2.204
            #
            # The snippet below prints weights and corresponding bond_factors.
            # [f"{w:.2f}: {1-np.log(w):.4f}" for w in np.linspace(0.3, 1, 25)]
            self.bond_factor = -math.log(self.min_weight) + 1
        self.log(f"Using a factor of {self.bond_factor:.6f} for bond detection.")

        # Set up primitive coordinate indices
        if typed_prims is None:
            self.set_primitive_indices(
                self.atoms,
                self.coords3d,
            )
        # Use supplied typed_prims
        else:
            self.log(f"{len(typed_prims)} primitives were supplied. Checking them.")
            valid_typed_prims = check_typed_prims(
                self.coords3d,
                typed_prims,
                bend_min_deg=self.bend_min_deg,
                dihed_max_deg=self.dihed_max_deg,
                lb_min_deg=self.lb_min_deg,
            )
            self.log(
                f"{len(valid_typed_prims)} primitives are valid at the current Cartesians."
            )
            self.typed_prims = valid_typed_prims
            self.set_inds_from_typed_prims(self.typed_prims)

        self.primitives = get_primitives(
            self.coords3d,
            self.typed_prims,
            logger=self.logger,
        )
        if self.bonds_only:
            self.bending_indices = list()
            self.dihedral_indices = list()
            self.linear_bend_indices = list()
            self.primitives = [
                prim for prim in self.primitives if isinstance(prim, Stretch)
            ]
        check_primitives(self.coords3d, self.primitives, logger=self.logger)

        self._prim_internals = self.eval(self.coords3d)
        self._prim_coords = np.array(
            [prim_int.val for prim_int in self._prim_internals]
        )

        bonds = len(self.bond_indices)
        bends = len(self.bending_indices) + len(self.linear_bend_indices)
        dihedrals = len(self.dihedral_indices)
        assert bonds + bends + dihedrals == len(self.primitives)
        self._bonds_slice = slice(bonds)
        self._bends_slice = slice(bonds, bonds + bends)
        self._dihedrals_slice = slice(bonds + bends, bonds + bends + dihedrals)
        self.backtransform_counter = 0

    def log(self, message):
        self.logger.debug(message)

    @property
    def coords3d(self):
        return self._coords3d

    @coords3d.setter
    def coords3d(self, coords3d):
        self._coords3d = coords3d.reshape(-1, 3)
        self._B_prim = None
        self._prim_coords = None
        self._prim_internals = None

    @property
    def primitives(self):
        return self._primitives

    @primitives.setter
    def primitives(self, primitives):
        self._primitives = primitives

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
        return np.array([prim_int.indices for prim_int in self.prim_internals[slice_]])

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

    def pinv(self, array, rcond=None):
        if rcond is None:
            rcond = self.rcond
        return np.linalg.pinv(array, rcond=rcond)

    @property
    def Bt_inv_prim(self):
        """Transposed generalized inverse of the primitive Wilson B-Matrix."""
        B = self.B_prim
        return self.pinv(B.dot(B.T)).dot(B)

    @property
    def Bt_inv(self):
        """Transposed generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return self.pinv(B.dot(B.T)).dot(B)

    @property
    def B_inv_prim(self):
        """Generalized inverse of the primitive Wilson B-Matrix."""
        B = self.B_prim
        return B.T.dot(self.pinv(B.dot(B.T)))

    @property
    def B_inv(self):
        """Generalized inverse of the Wilson B-Matrix."""
        B = self.B
        return B.T.dot(self.pinv(B.dot(B.T)))

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

        size_ = self.coords3d.size
        if int_gradient is None:
            return np.zeros((size_, size_))

        K_flat = np.zeros(size_ * size_)
        coords3d = self.coords3d
        for primitive, int_grad_item in zip(self.primitives, int_gradient):
            # Contract with gradient
            val = np.rad2deg(primitive.calculate(coords3d))
            self.log(f"K, {primitive}={val:.2f}°")
            # The generated code (d2q_d) seems unstable for these values...
            if isinstance(primitive, Torsion) and ((abs(val) < 1) or (abs(val) > 179)):
                self.log(f"Skipped 2nd derivative of {primitive} with val={val:.2f}°")
                continue
            # 2nd derivative of normal, but linear, bends is undefined.
            try:
                dg = int_grad_item * primitive.jacobian(coords3d)
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
            cart_inds = list(
                it.chain(*[range(3 * i, 3 * i + 3) for i in primitive.indices])
            )
            flat_inds = [
                row * size_ + col for row, col in it.product(cart_inds, cart_inds)
            ]
            K_flat[flat_inds] += dg
        K = K_flat.reshape(size_, size_)
        return K

    def log_int_grad_msg(self, int_gradient):
        if int_gradient is None:
            self.log(
                "Supplied 'int_gradient' is None. K matrix will be zero, "
                "so derivatives of the\nWilson-B-matrix are neglected in "
                "Hessian transformation."
            )

    def transform_hessian(self, cart_hessian, int_gradient=None):
        """Transform Cartesian Hessian to internal coordinates."""
        self.log_int_grad_msg(int_gradient)
        K = self.get_K_matrix(int_gradient)
        return self.Bt_inv_prim.dot(cart_hessian - K).dot(self.B_inv_prim)

    def backtransform_hessian(self, redund_hessian, int_gradient=None):
        """Transform Hessian in internal coordinates to Cartesians."""
        self.log_int_grad_msg(int_gradient)
        K = self.get_K_matrix(int_gradient)
        return self.B.T.dot(redund_hessian).dot(self.B) + K

    def project_hessian(self, H, shift=1000):
        """Expects a hessian in internal coordinates. See Eq. (11) in [1]."""
        P = self.P
        return P.dot(H).dot(P) + shift * (np.eye(P.shape[0]) - P)

    def project_vector(self, vector):
        """Project supplied vector onto range of B."""
        return self.P.dot(vector)

    def set_inds_from_typed_prims(self, typed_prims):
        linear_bend_types = (PrimTypes.LINEAR_BEND, PrimTypes.LINEAR_BEND_COMPLEMENT)
        per_type = {
            2: list(),
            3: list(),
            4: list(),
            "linear_bend": list(),
            "hydrogen_bond": list(),
        }
        for type_, *indices in typed_prims:
            key = len(indices)
            if type_ in (linear_bend_types):
                key = "linear_bend"
            per_type[key].append(indices)

            # Also keep hydrogen bonds
            if type_ == PrimTypes.HYDROGEN_BOND:
                per_type["hydrogen_bond"].append(indices)

        self.bond_indices = per_type[2]
        self.bending_indices = per_type[3]
        self.dihedral_indices = per_type[4]
        self.linear_bend_indices = per_type["linear_bend"]
        self.hydrogen_bond_indices = per_type["hydrogen_bond"]

        # TODO
        # self.fragments = coord_info.fragments

    def set_primitive_indices(
        self,
        atoms,
        coords3d,
    ):
        coord_info = setup_redundant(
            atoms,
            coords3d,
            factor=self.bond_factor,
            define_prims=self.define_prims,
            min_deg=self.bend_min_deg,
            dihed_max_deg=self.dihed_max_deg,
            lb_min_deg=self.lb_min_deg,
            min_weight=self.min_weight if self.weighted else None,
            logger=self.logger,
        )

        self.typed_prims = coord_info.typed_prims
        self.set_inds_from_typed_prims(self.typed_prims)

        self.fragments = coord_info.fragments

    def eval(self, coords3d, attr=None):
        prim_internals = eval_primitives(coords3d, self.primitives)

        if attr is not None:
            return np.array(
                [getattr(prim_internal, attr) for prim_internal in prim_internals]
            )

        return prim_internals

    def transform_int_step(self, int_step, pure=False):
        self.log(f"Backtransformation {self.backtransform_counter}")
        new_prim_internals, cart_step, failed = transform_int_step(
            int_step,
            self.coords3d.flatten(),
            self.prim_coords,
            self.B_prim,
            self.primitives,
            check_dihedrals=self.rebuild,
            rcond=self.rcond,
            logger=self.logger,
        )
        # Update coordinates
        if not pure:
            self.coords3d += cart_step.reshape(-1, 3)
            self.prim_internals = new_prim_internals
            self.backtransform_counter += 1
        return cart_step

    def __str__(self):
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"
