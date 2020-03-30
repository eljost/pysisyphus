import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Bend(Primitive):

    def __init__(self, *args, linear=False, complement=False, **kwargs):
        kwargs["calc_kwargs"] = ("complement", )
        super().__init__(*args, **kwargs)

        self.linear = linear
        self.complement = complement

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm

        # Eq. (24) in [1]
        uv_parallel = Bend.parallel(u, v)
        # Determine second vector for the cross product, to get an
        # orthogonal direction.
        if not uv_parallel:
            cross_vec = v
        elif uv_parallel and not Bend.parallel(u, (1, -1, 1)):
            cross_vec = (1, -1, 1)
        # elif not Bend.parallel(u, (1, 0, 0)) and (not Bend.parallel(v, (1, 0, 0))):
            # cross_vec = (1, 0, 0)
        else:
            cross_vec = (-1, 1, 1)
            # cross_vec = (0., 0., 1)
        w_dash = np.cross(u, cross_vec)
        w = w_dash / np.linalg.norm(w_dash)

        if complement:
            angle_rad = np.arccos(u.dot(w)) + np.arccos(w.dot(v))
        else:
            angle_rad = np.arccos(u.dot(v))

        if gradient:
            if complement:
                # Create second vector, orthongal to w and u
                w = np.cross(u, w)
            # print("pysis w", w)
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
