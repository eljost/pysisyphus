from math import sin

import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords.derivatives import d2q_a
from pysisyphus.linalg import cross3, norm3


class Bend(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        m, o, n = indices
        rho_mo = Bend.rho(atoms, coords3d, (m, o))
        rho_on = Bend.rho(atoms, coords3d, (o, n))
        rad = Bend._calculate(coords3d, indices)
        return (rho_mo * rho_on)**0.5 * (f_damping + (1-f_damping)*sin(rad))

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = norm3(u_dash)
        v_norm = norm3(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm

        udv = max(-1.0, min(1.0, u.dot(v)))
        angle_rad = np.arccos(udv)

        if gradient:
            cross_vec1 = ( 1, -1, 1)
            cross_vec2 = (-1,  1, 1)

            # Determine second vector for the cross product, to get an
            # orthogonal direction. Eq. (24) in [1]
            uv_parallel = Bend.parallel(u, v)
            if not uv_parallel:
                cross_vec = v
            elif not Bend.parallel(u, cross_vec1):
                cross_vec = cross_vec1
            else:
                cross_vec = cross_vec2

            w_dash = cross3(u, cross_vec)
            w = w_dash / norm3(w_dash)

            uxw = cross3(u, w)
            wxv = cross3(w, v)

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

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_a(*coords3d[indices].flatten())
