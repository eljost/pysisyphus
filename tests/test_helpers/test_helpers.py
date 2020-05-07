import numpy as np
from scipy.spatial.transform import Rotation

from pysisyphus.helpers import align_coords
from pysisyphus.helpers import geom_loader


def get_coords(translate=0., euler=None):
    geom = geom_loader("lib:benzene.xyz")
    ref = geom.coords3d.copy()

    geom_ = geom.copy()
    geom_.coords += translate
    mod = geom_.coords

    if euler:
        rot = Rotation.from_euler("XYZ", euler, degrees=True)
        rot_mat = rot.as_matrix()
        mod = np.array([c.dot(rot_mat) for c in mod.reshape(-1, 3)])

    all_coords = (ref, mod)
    return all_coords


def test_align_coords_trans():
    all_coords = get_coords(translate=5.)
    aligned = align_coords(all_coords)
    np.testing.assert_allclose(aligned[0], aligned[1])


def test_align_coords_trans_rot():
    all_coords = get_coords(translate=5., euler=(90, 45, 17))
    aligned = align_coords(all_coords)
    np.testing.assert_allclose(aligned[0], aligned[1])
