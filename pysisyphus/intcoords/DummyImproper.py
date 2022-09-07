import numpy as np

from pysisyphus.intcoords.Torsion import Torsion


class DummyImproper(Torsion):
    def __init__(self, indices, *args, fix_inner=True, **kwargs):
        self.fix_inner = fix_inner
        kwargs["calc_kwargs"] = ("fix_inner",)
        super().__init__(indices, *args, **kwargs)
        self.log("DummyImproper is never checked for collinear atoms!")

    @staticmethod
    def _get_dummy_coords(coords3d, indices, r=1.889):
        r"""
    dummy -> D\
    atomy    | \
             |  \
         A---B---C
        """
        a, b, c = indices
        # Center
        B = coords3d[b]
        # Bond, pointing away from a to c
        AC_ = coords3d[c] - coords3d[a]
        AC = AC_ / np.linalg.norm(AC_)

        w = DummyImproper._get_cross_vec(coords3d, indices).flatten()
        dummy_vec = (np.eye(3) - np.outer(AC, AC)) @ w
        dummy_coords = B + r * dummy_vec
        return dummy_coords

    @staticmethod
    def get_coords3d_and_indices_ext(coords3d, indices):
        dummy_coords = DummyImproper._get_dummy_coords(coords3d, indices)
        dummy_ind = len(coords3d)
        coords3d_ext = np.zeros((len(coords3d) + 1, 3))
        coords3d_ext[:dummy_ind] = coords3d
        coords3d_ext[dummy_ind] = dummy_coords
        a, b, c = indices
        indices_ext = a, b, dummy_ind, c
        return coords3d_ext, indices_ext

    @staticmethod
    def _weight(*args, **kwargs):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, fix_inner=False):
        coords3d_ext, indices_ext = DummyImproper.get_coords3d_and_indices_ext(
            coords3d, indices
        )

        results = Torsion._calculate(coords3d_ext, indices_ext, gradient=gradient)
        if gradient:
            val, grad = results

            # import pdb; pdb.set_trace()  # fmt: skip
            # Remove entries that belong to the dummy atom.
            grad = grad[:-3]  # .reshape[](-1, 3)[[0, 1, 3]]

            # Zero out contributions of the inner two atoms in the torsion.
            # So basically only the atom at the first index moves.
            #
            # This usually degrades optimization convergence.
            # if fix_inner:
            # grad.reshape(-1, 3)[indices_ext[1:3]] = 0.0
            return val, grad.flatten()
        else:
            return results
