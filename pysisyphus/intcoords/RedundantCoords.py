# [1] https://doi.org/10.1063/1.1515483 optimization review
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates
# [3] https://doi.org/10.1016/0009-2614(95)00646-L lindh model hessian
# [4] 10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
#     Handling of corner cases
# [5] https://doi.org/10.1063/1.462844 , Pulay 1992

import itertools as it
import math
from operator import itemgetter

import numpy as np

from pysisyphus.config import (
        BEND_MIN_DEG,
        LB_MIN_DEG,
        DIHED_MAX_DEG,
)
from pysisyphus.linalg import svd_inv
from pysisyphus.intcoords.exceptions import PrimitiveNotDefinedException
from pysisyphus.intcoords.update import transform_int_step
from pysisyphus.intcoords.eval import (
    eval_primitives,
    check_primitives,
)

from pysisyphus.intcoords.logging_conf import logger
from pysisyphus.intcoords.PrimTypes import (
    normalize_prim_inputs,
    PrimTypes,
    # PrimType classes
    Bonds,
    Bends,
    DummyCoords,
    LinearBends,
    Cartesians,
    Dihedrals,
    OutOfPlanes,
    Rotations,
    Translations,
)

from pysisyphus.intcoords.setup import (
    setup_redundant,
    get_primitives,
)
from pysisyphus.intcoords.valid import check_typed_prims


class RedundantCoords:
    def __init__(
        self,
        atoms,
        coords3d,
        masses=None,
        bond_factor=1.3,
        typed_prims=None,
        define_prims=None,
        constrain_prims=None,
        freeze_atoms=None,
        freeze_atoms_exclude=False,
        define_for=None,
        bonds_only=False,
        check_bends=True,
        rebuild=True,
        bend_min_deg=BEND_MIN_DEG,
        dihed_max_deg=DIHED_MAX_DEG,
        lb_min_deg=LB_MIN_DEG,
        weighted=False,
        min_weight=0.3,
        # Corresponds to a threshold of 1e-7 for eigenvalues of G, as proposed by
        # Pulay in [5].
        svd_inv_thresh=3.16e-4,
        recalc_B=False,
        tric=False,
        hybrid=False,
        hbond_angles=False,
    ):
        self.atoms = atoms
        self.coords3d = np.reshape(coords3d, (-1, 3)).copy()
        self.masses = masses
        self.bond_factor = bond_factor
        if typed_prims is not None:
            typed_prims = normalize_prim_inputs(typed_prims)
        # Define additional primitives
        if define_prims is None:
            define_prims = list()
        self.define_prims = normalize_prim_inputs(define_prims)
        if freeze_atoms is None:
            freeze_atoms = list()
        self.freeze_atoms = np.array(freeze_atoms, dtype=int)
        self.freeze_atoms_exclude = freeze_atoms_exclude
        self.define_for = define_for
        # Constrain primitives
        if constrain_prims is None:
            constrain_prims = list()
        self.constrain_prims = normalize_prim_inputs(constrain_prims)
        self.bonds_only = bonds_only
        self.check_bends = check_bends
        self.rebuild = rebuild
        self.bend_min_deg = bend_min_deg
        self.dihed_max_deg = dihed_max_deg
        self.lb_min_deg = lb_min_deg
        self.weighted = weighted
        self.min_weight = float(min_weight)
        assert self.min_weight > 0.0, "min_weight must be a positive rational!"
        self.svd_inv_thresh = svd_inv_thresh
        self.recalc_B = recalc_B
        self.tric = tric
        self.hybrid = hybrid
        self.hbond_angles = hbond_angles

        self._B_prim = None
        # Lists for the other types of primitives will be created afterwards.
        self.logger = logger

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
        self.log(f"Using svd_inv_thresh={self.svd_inv_thresh:.4e} for inversions.")

        # Set up primitive coordinate indices
        if typed_prims is None:
            self.set_primitive_indices(
                self.atoms,
                self.coords3d,
            )
        # Use supplied typed_prims
        else:
            self.typed_prims = typed_prims + self.define_prims

        if self.bonds_only:
            self.typed_prims = self.bond_typed_prims

        self.primitives = get_primitives(
            self.coords3d,
            self.typed_prims,
            logger=self.logger,
        )

        # First evaluation of internal coordinates
        self._prim_internals = self.eval(self.coords3d)
        self._prim_coords = np.array(
            [prim_int.val for prim_int in self._prim_internals]
        )
        check_primitives(
            self.coords3d, self.primitives, B=self.B_prim, logger=self.logger
        )

        ref_num = len(self.typed_prims)
        if self.bonds_only:
            ref_num = len(self.bond_indices)
        assert len(self.primitives) == ref_num

        self.backtransform_counter = 0

    def set_inds_from_typed_prims(self, typed_prims):
        # These lists will hold the index of the respective typed_prims
        # in 'self.typed_prims'.
        self._bond_inds = list()
        self._bend_inds = list()
        self._linear_bend_inds = list()
        self._dihedral_inds = list()
        self._rotation_inds = list()
        self._translation_inds = list()
        self._cartesian_inds = list()
        self._outofplane_inds = list()
        self._dummycoord_inds = list()
        self._cartesian_inds = list()

        self._bond_atom_inds = list()
        self._bend_atom_inds = list()
        self._dihedral_atom_inds = list()

        self._bond_typed_prims = list()

        for i, (pt, *indices) in enumerate(typed_prims):
            if pt in Bonds:
                append_to = self._bond_inds
                self._bond_atom_inds.append(indices)
                self._bond_typed_prims.append((pt, *indices))
            elif pt in Bends:
                append_to = self._bend_inds
                self._bend_atom_inds.append(indices)
            elif pt in LinearBends:
                append_to = self._linear_bend_inds
            elif pt in Dihedrals:
                append_to = self._dihedral_inds
                self._dihedral_atom_inds.append(indices)
            elif pt in Rotations:
                append_to = self._rotation_inds
            elif pt in Translations:
                append_to = self._translation_inds
            elif pt in Cartesians:
                append_to = self._cartesian_inds
            elif pt in OutOfPlanes:
                append_to = self._outofplane_inds
            elif pt in DummyCoords:
                append_to = self._dummycoord_inds
            elif pt in Cartesians:
                append_to = self._cartesian_inds
            else:
                raise Exception("Unhandled PrimType!")
            append_to.append(i)

    def log(self, message):
        self.logger.debug(message)

    def clear(self):
        self._B_prim = None
        self._prim_coords = None
        self._prim_internals = None
        self._P = None

    @property
    def coords3d(self):
        return self._coords3d

    @coords3d.setter
    def coords3d(self, coords3d):
        self._coords3d = coords3d.reshape(-1, 3)
        self.clear()

    @property
    def typed_prims(self):
        return self._typed_prims

    @typed_prims.setter
    def typed_prims(self, typed_prims):
        self.log(f"Checking {len(typed_prims)} supplied typed primitives.")
        valid_typed_prims = check_typed_prims(
            self.coords3d,
            typed_prims,
            bend_min_deg=self.bend_min_deg,
            dihed_max_deg=self.dihed_max_deg,
            lb_min_deg=self.lb_min_deg,
            check_bends=self.check_bends,
        )

        def tp_sort(tp):
            pt, *indices = tp
            key = pt
            # We use the fact that list.sort is stable, that is elements that compare
            # equal retain their order. So we assign PrimTypes.ROTATION to all rotations,
            # to remain the ABC-order for each fragment. The same goes for the translations.
            if pt in Rotations:
                key = PrimTypes.ROTATION
            elif pt in Translations:
                key = PrimTypes.TRANSLATION
            elif pt in Cartesians:
                key = PrimTypes.CARTESIAN
            return key

        # Sort by PrimType
        valid_typed_prims.sort(key=tp_sort)

        self.log(
            f"{len(valid_typed_prims)} primitives are valid at the current Cartesians."
        )
        if len(valid_typed_prims) != len(typed_prims):
            self.log("Invalid primitives:")
            for i, invalid_prim in enumerate(set(typed_prims) - set(valid_typed_prims)):
                self.log(f"\t{i:02d}: {invalid_prim}")
        self._typed_prims = valid_typed_prims
        self.set_inds_from_typed_prims(self.typed_prims)

    @property
    def primitives(self):
        return self._primitives

    @primitives.setter
    def primitives(self, primitives):
        self._primitives = primitives

    @property
    def prim_indices_set(self):
        return set([tuple(indices) for pt, *indices in self.typed_prims])

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

    def get_prim_internals_by_indices(self, indices):
        if len(indices) == 0:
            pis = []
        elif len(indices) == 1:
            pis = [self.prim_internals[indices[0]]]
        else:
            pis = itemgetter(*indices)(self.prim_internals)
        return pis

    @property
    def bond_indices(self):
        return self._bond_inds

    @property
    def bond_atom_indices(self):
        return self._bond_atom_inds

    @property
    def bond_typed_prims(self):
        return self._bond_typed_prims

    @property
    def bend_indices(self):
        return self._bend_inds

    @property
    def bend_atom_indices(self):
        return self._bend_atom_inds

    @property
    def linear_bend_indices(self):
        return self._linear_bend_inds

    @property
    def dihedral_indices(self):
        return self._dihedral_inds

    @property
    def dihedral_atom_indices(self):
        return self._dihedral_atom_inds

    @property
    def rotation_indices(self):
        return self._rotation_inds

    @property
    def translation_indices(self):
        return self._translation_inds

    @property
    def cartesian_indices(self):
        return self._cartesian_inds

    @property
    def outofplane_indices(self):
        return self._outofplane_inds

    @property
    def coords(self):
        return self.prim_coords

    def get_index_of_typed_prim(self, typed_prim):
        """Index in self.typed_prims for the supplied typed_prim."""
        ref_len = len(typed_prim)
        ref_inds = typed_prim[1:]
        for i, tp in enumerate(self.typed_prims):
            if (len(tp) != ref_len) or tp[0] != typed_prim[0]:
                continue

            if (tp[1:] == ref_inds) or (tp[1:] == ref_inds[::-1]):
                return i
        self.log(f"Typed primitive {typed_prim} is not defined!")
        raise PrimitiveNotDefinedException(typed_prim)

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

    def inv_B(self, B):
        return B.T.dot(svd_inv(B.dot(B.T), thresh=self.svd_inv_thresh, hermitian=True))

    def inv_Bt(self, B):
        return svd_inv(B.dot(B.T), thresh=self.svd_inv_thresh, hermitian=True).dot(B)

    @property
    def Bt_inv_prim(self):
        """Transposed generalized inverse of the primitive Wilson B-Matrix."""
        return self.inv_Bt(self.B_prim)

    @property
    def Bt_inv(self):
        """Transposed generalized inverse of the Wilson B-Matrix."""
        return self.inv_Bt(self.B)

    @property
    def B_inv_prim(self):
        """Generalized inverse of the primitive Wilson B-Matrix."""
        return self.inv_B(self.B_prim)

    @property
    def B_inv(self):
        """Generalized inverse of the Wilson B-Matrix."""
        return self.inv_B(self.B)

    @property
    def constrained_indices(self):
        return [self.typed_prims.index(cp) for cp in self.constrain_prims]

    @property
    def C(self):
        """Diagonal matrix. Entries for constraints are set to one."""
        size = len(self.typed_prims)
        C = np.zeros((size, size))
        inds = self.constrained_indices
        C[inds, inds] = 1
        return C

    @property
    def P(self):
        """Projection matrix onto B. See [1] Eq. (4)."""
        if self._P is None:
            P = self.B.dot(self.B_inv)
            # Modify projector, so constrained coordinates are projected out.
            if self.constrain_prims:
                C = self.C
                CPC_inv = svd_inv(C.dot(P).dot(C), thresh=self.svd_inv_thresh)
                P = P - P.dot(C).dot(CPC_inv).dot(C).dot(P)
            self._P = P
        return self._P

    def transform_forces(self, cart_forces):
        """Combination of Eq. (9) and (11) in [1]."""
        return self.P.dot(self.Bt_inv.dot(cart_forces))

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
            try:
                dg = int_grad_item * primitive.jacobian(coords3d)
            # 2nd derivative of normal, but linear, bends is undefined.
            except (ValueError, ZeroDivisionError):
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
            tric=self.tric,
            hybrid=self.hybrid,
            hbond_angles=self.hbond_angles,
            freeze_atoms=self.freeze_atoms if self.freeze_atoms_exclude else None,
            define_for=self.define_for,
            logger=self.logger,
        )

        self.typed_prims = coord_info.typed_prims
        for cp in self.constrain_prims:
            if cp not in self.typed_prims:
                self.typed_prims.append(cp)

        self.fragments = coord_info.fragments

    def eval(self, coords3d, attr=None):
        prim_internals = eval_primitives(coords3d, self.primitives)

        if attr is not None:
            return np.array(
                [getattr(prim_internal, attr) for prim_internal in prim_internals]
            )

        return prim_internals

    def transform_int_step(self, int_step, update_constraints=False, pure=False):
        self.log(f"Backtransformation {self.backtransform_counter}")

        def Bt_inv_prim_getter(cart_coords):
            coords3d = cart_coords.reshape(-1, 3)
            B_prim = np.zeros((len(self.primitives), coords3d.size))
            for i, primitive in enumerate(self.primitives):
                _, gradient = primitive.calculate(coords3d, gradient=True)
                B_prim[i] = gradient
            return self.inv_Bt(B_prim)

        new_prim_internals, cart_step, failed = transform_int_step(
            int_step,
            self.coords3d.flatten(),
            self.prim_coords,
            self.Bt_inv_prim,
            self.primitives,
            typed_prims=self.typed_prims,
            check_dihedrals=self.rebuild,
            check_bends=self.rebuild,
            bend_min_deg=self.bend_min_deg,
            bend_max_deg=self.lb_min_deg,
            freeze_atoms=self.freeze_atoms,
            constrained_inds=self.constrained_indices,
            update_constraints=update_constraints,
            logger=self.logger,
            Bt_inv_prim_getter=Bt_inv_prim_getter if self.recalc_B else None,
        )
        # Update coordinates
        if not pure:
            self.coords3d += cart_step.reshape(-1, 3)
            self.prim_internals = new_prim_internals
            self.backtransform_counter += 1
        return cart_step

    def print_typed_prims(self):
        for i, tp in enumerate(self.typed_prims):
            print(i, tp)

    def __str__(self):
        bonds = len(self.bond_indices)
        bends = len(self.bending_indices)
        dihedrals = len(self.dihedral_indices)
        name = self.__class__.__name__
        return f"{name}({bonds} bonds, {bends} bends, {dihedrals} dihedrals)"


class TRIC(RedundantCoords):
    def __init__(self, *args, **kwargs):
        kwargs["tric"] = True
        kwargs["recalc_B"] = True
        super().__init__(*args, **kwargs)


class HybridRedundantCoords(RedundantCoords):
    def __init__(self, *args, **kwargs):
        kwargs["hybrid"] = True
        super().__init__(*args, **kwargs)
