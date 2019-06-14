#!/usr/bin/env python3

import numpy as np

from pysisyphus.InternalCoordinates import RedundantCoords


class DLC(RedundantCoords):

    def set_delocalized_vectors(self, thresh=1e-6):
        """See [5] between Eq. (7) and Eq. (8) for advice regarding
        the threshold."""
        # Original Wilson B-matrix
        B = self.B
        G = B.dot(B.T)
        eigvals, eigvectors = np.linalg.eigh(G)
        #print(w)
        #print(w.shape)
        #print(v.T)
        nonzero_inds = np.abs(eigvals) > thresh
        sum_ = np.sum(nonzero_inds)
        active_eigvals = eigvals[nonzero_inds]
        active_vecs = eigvectors[:,nonzero_inds]
        # Transformation of B to the active coordinate set
        B = active_vecs.T.dot(B)
        Bt_inv = np.linalg.pinv(B.dot(B.T)).dot(B)
        # self.B = self.delocalized_vectors.T.dot(self.B_prim)
        # self.Bt_inv = np.linalg.pinv(self.B.dot(self.B.T)).dot(self.B)
        import pdb; pdb.set_trace()
        pass
