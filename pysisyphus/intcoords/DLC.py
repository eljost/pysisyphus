#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/10.1063/1.471864
#     The generation and use of delocalized internal coordinates in geometry
#     optimization.
#     Baker, 1996

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.linalg import gram_schmidt


class DLC(RedundantCoords):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.U = self.get_active_set(self.B_prim)
        # Keep a copy of the original active set, in case self.U gets
        # modified by applying some constraints.
        self.U_unconstrained = self.U.copy()

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, U):
        self._U = U
        # Needed for back-transformation to primitive internals
        self._Ut_inv = np.linalg.pinv(self.U.T)

    @property
    def Ut_inv(self):
        return self._Ut_inv

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

    def transform_int_step(self, step, *args, **kwargs):
        """As the transformation is done in primitive internal coordinates
        we convert the DLC back to primitive coordinates."""
        # Or: prim_step = (step*self.U).sum(axis=1)
        prim_step = self.Ut_inv.dot(step)
        return super().transform_int_step(prim_step, *args, **kwargs)

    def get_active_set(self, B, thresh=1e-6):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)

        nonzero_inds = np.abs(eigvals) > thresh
        active_eigvals = eigvals[nonzero_inds]
        return eigvectors[:,nonzero_inds]

    def project_primitive_on_active_set(self, prim_ind):
        prim_vec = np.zeros(self.U.shape[0])
        prim_vec[prim_ind] = 1
        c_proj = (np.einsum("i,ij->j", prim_vec, self.U) * self.U).sum(axis=1)
        c_proj /= np.linalg.norm(c_proj)
        return c_proj

    def apply_constraints(self, constraint_vecs):
        original_U_shape = self.U_unconstrained.shape
        constr_num = constraint_vecs.shape[1]
        V = np.concatenate((constraint_vecs, self.U), axis=1)
        orthonormalized = gram_schmidt(V.T).T
        # During Gram-Schmidt a number of columns of U should have
        # dropped out. They are replaced by the constraint_vecs, so
        # the total shape of should not change.
        assert orthonormalized.shape[1] == original_U_shape[1]
        # Remove constraint vectors
        # [1] states that we somehow have to keep the constraint vectors
        # (or the corresponding unit vectors) for the iterative
        # back-transformation. Right now I don't understand why we would
        # have to do this ([1], p. 10 (200), right column).
        U_proj = orthonormalized[:,constr_num:]
        return U_proj

    def freeze_primitives(self, prim_atom_indices):
        """Freeze primitive internal coordinates.

        Parameters
        ----------
        prim_atom_indices : iterable of atom index iterables
            Iterable containing atom index iterables that define the primitive
            internal to be frozen.
        """
        prim_indices = [self.get_index_of_prim_coord(pai)
                        for pai in prim_atom_indices
        ]
        _ = [self.project_primitive_on_active_set(pi)
                                for pi in prim_indices
        ]
        projected_primitives = np.array(_).T
        U_proj = self.apply_constraints(projected_primitives)
        self.U = U_proj
