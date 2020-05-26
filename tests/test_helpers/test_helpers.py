import numpy as np
from scipy.spatial.transform import Rotation

from pysisyphus.helpers import align_coords, geom_loader
from pysisyphus.interpolate import interpolate


def get_geoms(translate=0., euler=None):
    geom = geom_loader("lib:benzene.xyz")

    geom_mod = geom.copy()
    geom_mod.coords += translate
    mod = geom_mod.coords

    if euler:
        rot = Rotation.from_euler("XYZ", euler, degrees=True)
        rot_mat = rot.as_matrix()
        mod = np.array([c.dot(rot_mat) for c in mod.reshape(-1, 3)])
        geom_mod.coords = mod.flatten()

    return (geom, geom_mod)


def test_align_coords_trans():
    geoms = get_geoms(translate=5.)
    all_coords = [geom.coords for geom in geoms]
    aligned = align_coords(all_coords)
    np.testing.assert_allclose(aligned[0], aligned[1])


def test_align_coords_trans_rot():
    geoms = get_geoms(translate=5., euler=(90, 45, 17))
    all_coords = [geom.coords for geom in geoms]
    aligned = align_coords(all_coords)
    np.testing.assert_allclose(aligned[0], aligned[1])


def test_align_coords_trans_rot_3d():
    geoms = get_geoms(translate=5., euler=(90, 45, 17))
    all_coords = [geom.coords3d for geom in geoms]
    aligned = align_coords(all_coords)
    np.testing.assert_allclose(aligned[0], aligned[1])


def test_align_coords_interpolate():
    geoms = get_geoms(translate=5., euler=(0., 0., 90.))
    interpolated = interpolate(*geoms, 10, kind="lst")
    all_coords = [geom.coords for geom in interpolated]
    aligned = align_coords(all_coords)

    # from pysisyphus.xyzloader import coords_to_trj
    # trj_fn = "interpol.trj"
    # atoms = geoms[0].atoms
    # coords_to_trj("interpol.trj", atoms, all_coords)
    # trj_fn = "aligned.trj"
    # coords_to_trj("aligned.trj", atoms, aligned)

    np.testing.assert_allclose(aligned[0], aligned[-1], atol=1e-10)
