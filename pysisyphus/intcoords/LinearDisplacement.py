import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords.derivatives import dq_ld, d2q_ld


class LinearDisplacement(Primitive):
    def __init__(self, *args, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", "cross_vec")
        super().__init__(*args, **kwargs)

        self.complement = complement
        self.cross_vec = None

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        raise Exception("Not yet implemented!")

    def calculate(self, coords3d, indices=None, gradient=False):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().calculate(coords3d, indices, gradient)

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False, cross_vec=None):
        m, o, n = indices
        w_dash = coords3d[n] - coords3d[m]
        w = w_dash / np.linalg.norm(w_dash)

        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u = u_dash / np.linalg.norm(u_dash)
        v = v_dash / np.linalg.norm(v_dash)

        # Vector for cross product to determine first orthogonal direction
        if cross_vec is None:
            cross_vec = LinearDisplacement._get_cross_vec(coords3d, indices)

        if complement:
            cross_vec = np.cross(w, cross_vec)
        cross_vec /= np.linalg.norm(cross_vec)

        # Orthogonal direction
        y = np.cross(w, cross_vec)
        y /= np.linalg.norm(y)

        lin_disp = y.dot(u) + y.dot(v)

        if gradient:
            row = np.zeros_like(coords3d)
            row[indices] = dq_ld(*coords3d[indices].flatten(), *cross_vec).reshape(-1, 3)
            return lin_disp, row.flatten()

        return lin_disp

    def jacobian(self, coords3d, indices=None):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().jacobian(coords3d, indices)

    @staticmethod
    def _jacobian(coords3d, indices, complement=False, cross_vec=None):
        if cross_vec is None:
            cross_vec = LinearDisplacement._get_cross_vec(coords3d, indices)

        if complement:
            m, _, n = indices
            w_dash = coords3d[n] - coords3d[m]
            w = w_dash / np.linalg.norm(w_dash)
            cross_vec = np.cross(w, cross_vec)

        return d2q_ld(*coords3d[indices].flatten(), *cross_vec)

    def __str__(self):
        return (
            f"LinearDisplacement({tuple(self.indices)}, complement={self.complement})"
        )
