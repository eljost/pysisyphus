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

        self.set_active_set()

    @property
    def U(self):
        return self._U

    @U.setter
    def U(self, U):
        self._U = U
        # Needed for back-transformation to primitive internals
        # For now we use a pseudo-inverse instead of a regular inverse,
        # as some columns of U may be zero from constraints (resulting in
        # singular U matrix that cannot be readily inverted).
        self._Ut_inv = np.linalg.pinv(self.U.T)

    @property
    def Ut_inv(self):
        return self._Ut_inv

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
        not_defined = [prim_coord for prim_coord, prim_ind
                       in zip(prim_atom_indices, prim_indices)
                       if prim_ind is None
        ]
        assert None not in prim_indices, \
            f"Some primitive internals are not defined ({not_defined})!"
        _ = [self.project_primitive_on_active_set(pi)
                                for pi in prim_indices
        ]
        projected_primitives = np.array(_).T
        self.constraints = projected_primitives
