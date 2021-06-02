import pytest

from pysisyphus.xyzloader import write_geoms_to_trj
from pysisyphus.helpers import geom_loader
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.interpolate.LST import LST
from pysisyphus.interpolate.IDPP import IDPP
from pysisyphus.interpolate.Redund import Redund


def test_idpp():
    initial = geom_loader("lib:09_htransfer_product.xyz")
    final = geom_loader("lib:10_po_diss_product_xtbopt.xyz")

    geoms = (initial, final)
    idpp = IDPP(geoms, 18, align=True)
    geoms = idpp.interpolate_all()
    # idpp.all_geoms_to_trj("idpp_opt.trj")

    assert len(geoms) == 20


@pytest.mark.parametrize(
    "interpol_cls",
    [
        Interpolator,
        LST,
        IDPP,
        Redund,
    ]
)
def test_ala_dipeptide_interpol(interpol_cls):
    initial = geom_loader("lib:dipeptide_init.xyz")
    final = geom_loader("lib:dipeptide_fin.xyz")

    geoms = (initial, final)
    interpolator = interpol_cls(geoms, 28, align=True)
    geoms = interpolator.interpolate_all()
    # interpolator.all_geoms_to_trj("interpolated.trj")

    assert len(geoms) == 30
