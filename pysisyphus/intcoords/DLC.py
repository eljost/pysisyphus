#!/usr/bin/env python3

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords


class DLC(RedundantCoords):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.set_active_set()

        self._Ut_inv = np.linalg.pinv(self.active_set.T)

    @property
    def U(self):
        return self.active_set

    @property
    def Ut_inv(self):
        return self._Ut_inv

    @property
    def coords(self):
        import pdb; pdb.set_trace()
        coords_org = super().coords
        return self.active_set.T.dot(coords_org)

    @property
    def B(self):
        """Wilson B-Matrix"""
        B_org = super().B
        print("Call DLC B")
        return self.active_set.T.dot(B_org)

    def project_hessian(self, H):
        """As we work in the non-redundant subspace we don't have to project
        the hessian."""
        return H

    def update_internals(self, new_cartesians, prev_internals):
        # Transform prev_internals to primitive internals
        prev_internals_prim = self.Ut_inv.dot(prev_internals)
        updated_internals = super().update_internals(new_cartesians, prev_internals_prim)
        return self.U.T.dot(updated_internals)

    def set_active_set(self, thresh=1e-6):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        # Original Wilson B-matrix
        B = super().B
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)

        nonzero_inds = np.abs(eigvals) > thresh
        # sum_ = np.sum(nonzero_inds)
        active_eigvals = eigvals[nonzero_inds]
        self.active_set = eigvectors[:,nonzero_inds]
        # B = active_vecs.T.dot(B)
        # Bt_inv = np.linalg.pinv(B.dot(B.T)).dot(B)
        # self.B = self.delocalized_vectors.T.dot(self.B_prim)
        # self.Bt_inv = np.linalg.pinv(self.B.dot(self.B.T)).dot(self.B)
