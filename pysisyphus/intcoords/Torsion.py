import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Torsion(Primitive):

    def __init__(self, *args, cos_tol=1e-9, **kwargs):
        kwargs["calc_kwargs"] = ("cos_tol", )
        super().__init__(*args, **kwargs)

        self.cos_tol = cos_tol

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, cos_tol=1e-9):
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

        # Restrict cos_dihed to [-1, 1]
        if cos_dihed >= 1 - cos_tol:
            dihedral_rad = 0
        elif cos_dihed <= -1 + cos_tol:
            dihedral_rad = np.arccos(-1)
        else:
            dihedral_rad = np.arccos(cos_dihed)

        if dihedral_rad != np.pi:
            # wxv = np.cross(w, v)
            # if wxv.dot(u) < 0:
            if vxw.dot(u) < 0:
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
