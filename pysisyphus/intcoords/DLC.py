#!/usr/bin/env python3

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.linalg import gram_schmidt


class DLC(RedundantCoords):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # U
        self._U = self.get_active_set(self.B_prim)
        # Needed for back-transformation to primitive internals
        self._Ut_inv = np.linalg.pinv(self.U.T)

    @property
    def U(self):
        return self._U

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
        """As we work in the non-redundant subspace we don't have to project
        the hessian."""
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

    def constrain_active_set(self, constraint_vecs):
        pass
