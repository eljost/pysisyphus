import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class LinearBend(Primitive):

    def __init__(self, *args, linear=False, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", )
        super().__init__(*args, **kwargs)

        self.linear = linear
        self.complement = complement

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False):
        m, o, n, p = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm

        cross_vec = coords3d[p] - coords3d[o]
        w_dash = np.cross(u, cross_vec)
        if complement:
            # import pdb; pdb.set_trace()
            w_dash = -np.cross(u, w_dash)
            print("uw comp", u.dot(w_dash))
            print("vw comp", v.dot(w_dash))
        w = w_dash / np.linalg.norm(w_dash)

        angle_rad = np.arccos(u.dot(w)) + np.arccos(w.dot(v))

        if gradient:
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

    def __str__(self):
        return f"Bend({tuple(self.indices)}, linear={self.linear}, " \
               f"complement={self.complement})"
