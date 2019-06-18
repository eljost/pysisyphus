#!/usr/bin/env python3

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords


class DLC(RedundantCoords):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # U
        self._U = self.set_active_set(self.B_prim)
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
        prim_step = self.Ut_inv.dot(step)
        return super().transform_int_step(prim_step, *args, **kwargs)

    def set_active_set(self, B, thresh=1e-6):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)

        nonzero_inds = np.abs(eigvals) > thresh
        active_eigvals = eigvals[nonzero_inds]
        return eigvectors[:,nonzero_inds]
