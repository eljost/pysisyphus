# [1] https://aip.scitation.org/doi/10.1063/1.471864
#     The generation and use of delocalized internal coordinates in geometry
#     optimization.
#     Baker, 1996

import numpy as np

from pysisyphus.intcoords import RedundantCoords
from pysisyphus.linalg import gram_schmidt


class DLC(RedundantCoords):
    def __init__(self, *args, full_set=True, **kwargs):
        super().__init__(*args, **kwargs)

        self.full_set = full_set
        self.set_active_set()

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, U):
        self._U = U

    @property
    def constraints(self):
        return self._constraints

    @constraints.setter
    def constraints(self, constraints):
        self.U = self.U_unconstrained.copy()
        constraints = np.array(constraints)
        constraints.flags.writeable = False
        # Constraint columns should be same length as columns of U.
        assert constraints.shape[0] == self.U.shape[0]
        self._constraints = constraints
        U_constrained = self.get_constrained_U(self._constraints)

        # # Replace the constraint-columns with zeros vectors, so the total
        # # number of coords doesn't change. This also zeros any force components
        # # that belong to the constraints.
        # zero_arr = np.zeros_like(constraints)
        # self.U = np.concatenate((zero_arr, U_constrained), axis=1)
        U_constrained = np.concatenate((constraints, U_constrained), axis=1)
        self.U = U_constrained

    def reset_constraints(self):
        self._constraints = tuple()
        self.U = self.U_unconstrained

    @property
    def coords(self):
        return self.U.T.dot(self.prim_coords)

    @property
    def B(self):
        """Wilson B-Matrix in the non-redundant subspace."""
        return self.U.T.dot(self.B_prim)

    def project_hessian(self, H):
        """As we already work in the non-redundant subspace we don't have
        to project/shift the hessian as we do it for simple redundant
        internal coordinates."""
        return H

    def transform_hessian(self, cart_hessian, int_gradient=None):
        """Transform Cartesian Hessian to DLC."""
        # Transform the DLC gradient to primitive coordinates
        if int_gradient is not None:
            prim_gradient = (self.U * int_gradient).sum(axis=1)
        else:
            prim_gradient = None
        H = super().transform_hessian(cart_hessian, prim_gradient)
        return self.U.T.dot(H).dot(self.U)

    def backtransform_hessian(self, *args, **kwargs):
        raise Exception("Check if we can just use the parents method.")

    def transform_int_step(self, step, *args, **kwargs):
        """As the transformation is done in primitive internal coordinates
        we convert the DLC back to primitive coordinates."""
        prim_step = (step * self.U).sum(axis=1)
        return super().transform_int_step(prim_step, *args, **kwargs)

    def get_active_set(self, B, inv_thresh=None):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        if self.weighted:
            weights = np.array(
                [prim.weight(self.atoms, self.coords3d) for prim in self.primitives]
            )
            self.log(
                "Weighting B-matrix:\n"
                f"\tWeights: {np.array2string(weights, precision=4)}\n"
                f"\tmax(weights)={weights.max():.4f}, "
                f"min(weights)={weights.min():.4f}, ({len(weights)} primitives)"
            )
            B = np.diag(weights).dot(B)

        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)

        if inv_thresh is None:
            # The singular values of G=B B^T correspond to the square roots of the
            # eigenvalues of G.
            #
            # S = sqrt(w)
            # w = S**2
            #
            # To stay consistent with the SVD, we derive the eigenvalue threshold from
            # the SVD threshold.
            inv_thresh = self.svd_inv_thresh ** 2

        if self.full_set:
            use_inds = np.full_like(eigvals, False, dtype=bool)
            dof = 3 * len(self.atoms) - 6
            use_inds[-dof:] = True
        else:
            use_inds = np.abs(eigvals) > inv_thresh
        use_eigvals = eigvals[use_inds]
        min_eigval = use_eigvals.min()
        self.log(
            f"Diagonalizing G yielded {use_inds.sum()} DLCs. Smallest eigenvalue "
            f"is {min_eigval:.4e}"
        )
        return eigvectors[:, use_inds]

    def set_active_set(self):
        self.U = self.get_active_set(self.B_prim)
        # Keep a copy of the original active set, in case self.U gets
        # modified by constraint application.
        self.U_unconstrained = self.U.copy()
        self.original_U_shape = self.U_unconstrained.shape
        self._constraints = tuple()

    def project_primitive_on_active_set(self, prim_ind):
        prim_vec = np.zeros(self.U.shape[0])
        prim_vec[prim_ind] = 1
        c_proj = (np.einsum("i,ij->j", prim_vec, self.U) * self.U).sum(axis=1)
        c_proj /= np.linalg.norm(c_proj)
        return c_proj

    def get_constrained_U(self, constraint_vecs):
        # Constraints are organized in columns
        constr_num = constraint_vecs.shape[1]
        V = np.concatenate((constraint_vecs, self.U), axis=1)
        orthonormalized = gram_schmidt(V.T).T
        # During Gram-Schmidt a number of columns of U should have
        # dropped out. They are replaced by the constraint_vecs, so
        # the total shape of should not change.
        assert orthonormalized.shape[1] == self.original_U_shape[1]
        # Remove constraint vectors
        # [1] states that we somehow have to keep the constraint vectors
        # (or the corresponding unit vectors) for the iterative
        # back-transformation. Right now I don't understand why we would
        # have to do this ([1], p. 10 (200), right column).
        U_proj = orthonormalized[:, constr_num:]
        return U_proj

    def freeze_primitives(self, typed_prims):
        """Freeze primitive internal coordinates.

        Parameters
        ----------
        typed_prims : iterable of typed primitives
            Iterable containing typed_primitives, starting with a PrimType and
            followed by atom indices.
        """
        prim_indices = [self.get_index_of_typed_prim(tp) for tp in typed_prims]
        not_defined = [
            tp for tp, prim_ind in zip(typed_prims, prim_indices) if prim_ind is None
        ]
        assert (
            None not in prim_indices
        ), f"Some primitive internals are not defined ({not_defined})!"
        projected_primitives = np.array(
            [self.project_primitive_on_active_set(pi) for pi in prim_indices]
        ).T
        self.constraints = projected_primitives


class HDLC(DLC):
    def __init__(self, *args, **kwargs):
        kwargs["hybrid"] = True
        super().__init__(*args, **kwargs)
