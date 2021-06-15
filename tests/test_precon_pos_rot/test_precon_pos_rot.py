import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.drivers.precon_pos_rot import precon_pos_rot
from pysisyphus.helpers import geom_loader


Bm = Benchmark("precon_pos_rot")


def test_precon_pos_rot_figure2(this_dir):
    educt, product = geom_loader("figure2_mod.trj")
    rgeom, pgeom = precon_pos_rot(educt, product, prefix="figure2_")


@pytest.mark.skip
@pytest.mark.parametrize("fn, geoms, ref_energy", Bm)
def test_birkholz_benchmark(fn, geoms, ref_energy):
    reactants, _, products = geoms
    rgeom, pgeom = precon_pos_rot(reactants, products)
