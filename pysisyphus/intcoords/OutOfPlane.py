import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords.derivatives import dq_oop, d2q_oop


class OutOfPlane(Primitive):
    """
    [1] https://doi.org/10.1002/(SICI)1096-987X(19990730)20:10<1067::AID-JCC9>3.0.CO;2-V
        Lee, 1999
    """

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        raise Exception("Not yet implemented!")

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        """
              P
            / | \
           /  |  \
          u'  v'  w'
         /    |    \
        m     n     o
        """
        # p is apex
        m, n, o, p = indices

        u_dash = coords3d[m] - coords3d[p]
        v_dash = coords3d[n] - coords3d[p]
        w_dash = coords3d[o] - coords3d[p]

        u = u_dash / np.linalg.norm(u_dash)
        v = v_dash / np.linalg.norm(v_dash)
        w = w_dash / np.linalg.norm(w_dash)

        z_dash = np.cross(u, v) + np.cross(v, w) + np.cross(w, u)
        z = z_dash / np.linalg.norm(z_dash)

        oop_coord = z.dot(u)

        if gradient:
            grad = dq_oop(*coords3d[m], *coords3d[n], *coords3d[o], *coords3d[p])
            grad = grad.reshape(4, 3)
            row = np.zeros_like(coords3d)
            row[m, :] = grad[0]
            row[n, :] = grad[1]
            row[o, :] = grad[2]
            row[p, :] = grad[3]
            row = row.flatten()
            return oop_coord, row

        return oop_coord

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_oop(*coords3d[indices].flatten())
