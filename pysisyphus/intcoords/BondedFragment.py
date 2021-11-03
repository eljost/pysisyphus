import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.linalg import norm3


class BondedFragment(Primitive):
    def __init__(self, indices, bond_indices, **kwargs):
        from_frag, to_ = bond_indices
        assert from_frag in indices
        assert to_ not in indices
        self.bond_indices = list(bond_indices)
        kwargs["calc_kwargs"] = ("bond_indices",)
        super().__init__(indices, **kwargs)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, bond_indices=None):
        from_frag, to_ = bond_indices
        bond = coords3d[to_] - coords3d[from_frag]
        value = norm3(bond)

        if gradient:
            bond_normed = bond / value
            row = np.zeros_like(coords3d)
            row[indices, :] = -bond_normed
            row = row.flatten()
            return value, row
        return value
