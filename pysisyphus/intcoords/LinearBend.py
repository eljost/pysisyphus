import numpy as np

from pysisyphus.intcoords.derivatives import dq_lb
from pysisyphus.intcoords.Primitive import Primitive

# [1] 10.1080/00268977200102361
#     Hoy, 1972
# [2] 10.1063/1.474377
#     Chuang, 1997
# [3] 10.1063/1.468630
#     Jackels, 1995
# Refs. [2] and [3] give a short discussion of the linear bends.


class LinearBend(Primitive):

    def __init__(self, *args, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", )
        super().__init__(*args, **kwargs)

        self.complement = complement

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        # v = v_dash / v_norm  # Never used ...

        # Select initial vector for cross product, similar to
        # geomeTRIC. It must NOT be parallel to u and/or v.
        x_dash = coords3d[n] - coords3d[m]
        x_norm = np.linalg.norm(x_dash)
        x = x_dash / x_norm
        cross_vecs = np.eye(3)
        min_ind = np.argmin([np.dot(cv, x)**2 for cv in cross_vecs])
        cross_vec = cross_vecs[min_ind]
        # Generate first orthogonal direction
        w_dash = np.cross(u, cross_vec)
        w = w_dash / np.linalg.norm(w_dash)

        # Generate second orthogonal direction
        if complement:
            w = np.cross(u, w)

        lb_rad = w.dot(np.cross(u_dash, v_dash)) / (u_norm * v_norm)

        if gradient:
            # Foruth argument is the orthogonal direction
            grad = dq_lb(*coords3d[m], *coords3d[o], *coords3d[n], *w)
            grad = grad.reshape(-1, 3)

            row = np.zeros_like(coords3d)
            row[m, :] = grad[0]
            row[o, :] = grad[1]
            row[n, :] = grad[2]
            row = row.flatten()
            return lb_rad, row
        return lb_rad

    def __str__(self):
        return f"LinearBend({tuple(self.indices)}, complement={self.complement})"
