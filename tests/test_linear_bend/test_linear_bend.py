import numpy as np

from pysisyphus.helpers import geom_loader


def calc_linear_bend(coords3d, angle_ind, grad=False):
    m, o, n = angle_ind
    u_dash = coords3d[m] - coords3d[o]
    v_dash = coords3d[n] - coords3d[o]
    u_norm = np.linalg.norm(u_dash)
    v_norm = np.linalg.norm(v_dash)
    u = u_dash / u_norm
    v = v_dash / v_norm

    def parallel(u, v, thresh=1e-4):
        dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return (1 - abs(dot)) < thresh

    uv_parallel = parallel(u, v)
    if not uv_parallel:
        cross_vec = v
    elif uv_parallel and not parallel(u, (1, -1, 1)):
        cross_vec = (1, -1, 1)
    else:
        cross_vec = (-1, 1, 1)

    # Define plane in which bending is measured through vector perpendicular to it
    #   e_A in [6]
    w_dash = np.cross(u, cross_vec)
    w = w_dash / np.linalg.norm(w_dash)
    # Define second vector, perpendicular to e_A and u
    #   e_B in [6]
    x_dash = np.cross(w, u)
    x = x_dash / np.linalg.norm(x_dash)

    lBA = np.arccos(u.dot(w)) + np.arccos(w.dot(v))
    import pdb; pdb.set_trace()
    return

    angle_rad = np.arccos(u.dot(v))
    if grad:
        # Eq. (24) in [1]
        if self.are_parallel(u, v, angle_ind):
            tmp_vec = np.array((1, -1, 1))
            par = self.are_parallel(u, tmp_vec) and self.are_parallel(v, tmp_vec)
            tmp_vec = np.array((-1, 1, 1)) if par else tmp_vec
            w_dash = np.cross(u, tmp_vec)
        else:
            w_dash = np.cross(u, v)
        w_norm = np.linalg.norm(w_dash)
        w = w_dash / w_norm
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


def test_hcn_linear_bend():
    geom = geom_loader("lib:hcn.xyz")
    # print(geom)
    # geom.jmol()
    # 179.5
    lb = (1, 0, 2)
    calc_linear_bend(geom.coords3d, lb)
