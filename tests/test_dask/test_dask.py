from distributed import LocalCluster

from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader

from pysisyphus.calculators import ORCA
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.testing import using


@using("orca")
def test_dask():
    with LocalCluster(n_workers=2, threads_per_worker=1) as cluster:
        address = cluster.scheduler_address

        geoms = geom_loader("lib:ala_dipeptide_iso_b3lyp_631gd_10_images.trj")

        for i, geom in enumerate(geoms):
            calc = ORCA("hf-3c", pal=1, calc_number=i)
            geom.set_calculator(calc)

        neb = NEB(geoms, scheduler=address)

        max_cycles = 1
        opt = SteepestDescent(neb, align=True, max_cycles=max_cycles)
        opt.run()

    assert opt.cur_cycle == (max_cycles - 1)
