from distributed import LocalCluster
import pytest

from pysisyphus.cos.NEB import NEB
from pysisyphus.drivers import run_opt
from pysisyphus.helpers import geom_loader

from pysisyphus.calculators import ORCA
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.testing import using


@pytest.mark.skip_ci
@using("orca")
def test_dask():
    with LocalCluster(n_workers=5, threads_per_worker=1) as cluster:
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


@pytest.mark.skip_ci
@pytest.mark.parametrize("is_external", (True, False))
@using("orca")
def test_external(is_external):
    n_workers = 5
    cos_kwargs = {}
    if is_external:
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
        cos_kwargs["scheduler"] = cluster.scheduler_address
    else:
        cos_kwargs["cluster"] = True
        cos_kwargs["cluster_kwargs"] = {
            "n_workers": 5,
        }

    geoms = geom_loader("lib:ala_dipeptide_iso_b3lyp_631gd_10_images.trj")

    calc_number = 0

    def calc_getter():
        nonlocal calc_number
        calc = ORCA("hf-3c", pal=1, calc_number=calc_number)
        calc_number += 1
        return calc

    neb = NEB(geoms, **cos_kwargs)
    ref_address = neb.scheduler
    print(f"{is_external=}, {ref_address=}, {neb._external_scheduler}")

    opt_result = run_opt(
        neb,
        calc_getter,
        "sd",
        {
            "max_cycles": 1,
        },
    )
    opt = opt_result.opt

    if is_external:
        cluster.close()
    assert neb.scheduler == ref_address
