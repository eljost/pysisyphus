import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Bend(Primitive):

    def __init__(self, indices, complement=False):
        super().__init__(indices)

        self.complement = complement

    def _calculate(self, coords3d, indices, gradient):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm

        def parallel(u, v, thresh=1e-4):
            dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
            return (1 - abs(dot)) < thresh

        # Eq. (24) in [1]
        uv_parallel = parallel(u, v)
        if not uv_parallel:
            cross_vec = v
        elif uv_parallel and not parallel(u, (1, -1, 1)):
            cross_vec = (1, -1, 1)
        else:
            cross_vec = (-1, 1, 1)
        w_dash = np.cross(u, cross_vec)
        w = w_dash / np.linalg.norm(w_dash)

        if self.complement:
            angle_rad = np.arccos(u.dot(w)) + np.arccos(w.dot(v))
        else:
            angle_rad = np.arccos(u.dot(v))

        if gradient:
            if self.complement:
                w = -(u + v) / np.linalg.norm(u + v)
            uxw = np.cross(u, w)
            wxv = np.cross(w, v)

            row = np.zeros_like(coords3d)
            #                  |  m  |  n  |  o  |
            # -----------------------------------
            # sign_factor(amo) |  1  |  0  | -1  | first_term
            # sign_factor(ano) |  0  |  1  | -1  | second_term
            first_term = uxw / u_norm
            second_term = wxv / v_norm
            row[m,:] = first_term
            row[o,:] = -first_term - second_term
            row[n,:] = second_term
            row = row.flatten()
            return angle_rad, row
        return angle_rad
