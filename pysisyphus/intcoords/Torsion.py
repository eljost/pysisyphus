from math import sin

import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords import Bend
from pysisyphus.intcoords.derivatives import d2q_d



class Torsion(Primitive):

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        m, o, p, n = indices
        rho_mo = Torsion.rho(atoms, coords3d, (m, o))
        rho_op = Torsion.rho(atoms, coords3d, (o, p))
        rho_pn = Torsion.rho(atoms, coords3d, (p, n))
        rad_mop = Bend._calculate(coords3d, (m, o, p))
        rad_opn = Bend._calculate(coords3d, (o, p, n))
        return (
            (rho_mo * rho_op * rho_pn)**(1/3)
            * (f_damping + (1-f_damping)*sin(rad_mop))
            * (f_damping + (1-f_damping)*sin(rad_opn))
        )

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        m, o, p, n = indices
        u_dash = coords3d[m] - coords3d[o]
        v_dash = coords3d[n] - coords3d[p]
        w_dash = coords3d[p] - coords3d[o]
        u_norm = np.linalg.norm(u_dash)
        v_norm = np.linalg.norm(v_dash)
        w_norm = np.linalg.norm(w_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        w = w_dash / w_norm
        phi_u = np.arccos(u.dot(w))
        phi_v = np.arccos(-w.dot(v))
        uxw = np.cross(u, w)
        vxw = np.cross(v, w)
        cos_dihed = uxw.dot(vxw)/(np.sin(phi_u)*np.sin(phi_v))
        # Restrict cos_dihed to the allowed interval for arccos [-1, 1]
        cos_dihed = min(1, max(cos_dihed, -1))

        dihedral_rad = np.arccos(cos_dihed)

        # Arccos only returns values between 0 and π, but dihedrals can
        # also be negative. This is corrected now.
        #
        # (v ⨯ w) · u will be < 0 when both vectors point in different directions.
        #
        #  M  --->   N
        #  ^        ^
        #   \      /
        #    u    v    positive dihedral, M rotates into N clockwise
        #     \  /     (v ⨯ w) · u > 0, keep positive sign
        #      OwP
        #              w points downward, into the screen plane.
        #              The vector resulting from the cross-product is easily
        #              visualized with your right hand.
        #
        #  M
        #   \
        #  | u
        #  |  \
        #  |   OwP     negative dihedral, M rotates into N counter-clockwise
        #  v  /        (v ⨯ w) · u < 0, invert dihedral sign
        #    v
        #   /
        #  N
        #
        if (dihedral_rad != np.pi) and (vxw.dot(u) < 0):
            dihedral_rad *= -1

        if gradient:
            row = np.zeros_like(coords3d)
            #                  |  m  |  n  |  o  |  p  |
            # ------------------------------------------
            # sign_factor(amo) |  1  |  0  | -1  |  0  | 1st term
            # sign_factor(apn) |  0  | -1  |  0  |  1  | 2nd term
            # sign_factor(aop) |  0  |  0  |  1  | -1  | 3rd term
            # sign_factor(apo) |  0  |  0  | -1  |  1  | 4th term
            sin2_u = np.sin(phi_u)**2
            sin2_v = np.sin(phi_v)**2
            first_term  = uxw/(u_norm*sin2_u)
            second_term = vxw/(v_norm*sin2_v)
            third_term  = uxw*np.cos(phi_u)/(w_norm*sin2_u)
            fourth_term = -vxw*np.cos(phi_v)/(w_norm*sin2_v)
            row[m,:] = first_term
            row[n,:] = -second_term
            row[o,:] = -first_term + third_term - fourth_term
            row[p,:] = second_term - third_term + fourth_term
            row = row.flatten()
            return dihedral_rad, row
        return dihedral_rad

    @staticmethod
    def _jacobian(coords3d, indices):
        sign = np.sign(Torsion._calculate(coords3d, indices))
        return sign * d2q_d(*coords3d[indices].flatten())
