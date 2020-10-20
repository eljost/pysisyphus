from math import sin

import numpy as np

from pysisyphus.intcoords.derivatives import dq_lb, d2q_lb
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
        kwargs["calc_kwargs"] = ("complement", "cross_vec")
        super().__init__(*args, **kwargs)

        self.complement = complement
        self.cross_vec = None

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        m, o, n = indices
        rho_mo = LinearBend.rho(atoms, coords3d, (m, o))
        rho_on = LinearBend.rho(atoms, coords3d, (o, n))

        # Repeated code to avoid import of intcoords.Bend
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        rad = np.arccos(u.dot(v))

        return (rho_mo * rho_on) ** 0.5 * (f_damping + (1 - f_damping) * sin(rad))

    @staticmethod
    def _get_orthogonal_direction(coords3d, indices, complement=False, cross_vec=None):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        u = u_dash / u_norm

        if cross_vec is None:
            cross_vec = LinearBend._get_cross_vec(coords3d, indices)
        # Generate first orthogonal direction
        w_dash = np.cross(u, cross_vec)
        w = w_dash / np.linalg.norm(w_dash)

        # Generate second orthogonal direction
        if complement:
            w = np.cross(u, w)
        return w

    def calculate(self, coords3d, indices=None, gradient=False):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().calculate(coords3d, indices, gradient)

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, complement=False, cross_vec=None):
        m, o, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w = LinearBend._get_orthogonal_direction(
            coords3d, indices, complement, cross_vec
        )

        lb_rad = w.dot(np.cross(u_dash, v_dash)) / (u_norm * v_norm)

        if gradient:
            # Fourth argument is the orthogonal direction
            row = np.zeros_like(coords3d)
            row[indices] = dq_lb(*coords3d[indices].flatten(), *w).reshape(-1, 3)
            return lb_rad, row.flatten()
        return lb_rad

    def jacobian(self, coords3d, indices=None):
        if self.cross_vec is None:
            self.set_cross_vec(coords3d, indices)

        return super().jacobian(coords3d, indices)

    @staticmethod
    def _jacobian(coords3d, indices, complement=False, cross_vec=None):
        if cross_vec is None:
            cross_vec = LinearBend._get_cross_vec(coords3d, indices)

        w = LinearBend._get_orthogonal_direction(
            coords3d, indices, complement, cross_vec
        )
        return d2q_lb(*coords3d[indices].flatten(), *w)

    def __str__(self):
        return f"LinearBend({tuple(self.indices)}, complement={self.complement})"
