import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Stretch(Primitive):

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        n, m = indices
        bond = coords3d[m] - coords3d[n]
        bond_length = np.linalg.norm(bond)
        if gradient:
            bond_normed = bond / bond_length
            row = np.zeros_like(coords3d)
            # 1 / -1 correspond to the sign factor [1] Eq. 18
            row[m,:] =  bond_normed
            row[n,:] = -bond_normed
            row = row.flatten()
            return bond_length, row
        return bond_length
