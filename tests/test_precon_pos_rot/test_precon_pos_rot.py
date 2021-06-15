from pysisyphus.helpers import geom_loader, align_coords
from pysisyphus.drivers.precon_pos_rot import precon_pos_rot


def test_precon_pos_rot_figure2(this_dir):
    educt, product = geom_loader("figure2_mod.trj")
    rgeom, pgeom = precon_pos_rot(educt, product)
