import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.drivers import run_perf, print_perf_results
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@using("pyscf")
@pytest.mark.parametrize(
    "pals, mems, pal_range, mem_range, repeat",
    [
        (None, 500, (1, 3, 1), None, 2),
        ((1, 2, 3), 500, None, None, 1),
    ],
)
def test_run_perf(pals, mems, pal_range, mem_range, repeat):
    geom = geom_loader("lib:h2o.xyz")

    calc_number = 0

    def calc_getter():
        nonlocal calc_number
        calc = PySCF(basis="321g", verbose=0, calc_number=calc_number)
        calc_number += 1
        return calc

    p_res = run_perf(
        geom, calc_getter, pals=pals, mems=mems, pal_range=pal_range, repeat=repeat
    )
    assert all([len(p_res[k]) == repeat for k in p_res.keys()])
    print_perf_results(p_res)


@using("pyscf")
@pytest.mark.parametrize("kind", ["energy", "forces", "hessian"])
def test_run_perf_kinds(kind):
    geom = geom_loader("lib:h2o.xyz")

    calc_number = 0

    def calc_getter():
        nonlocal calc_number
        calc = PySCF(basis="321g", verbose=0, calc_number=calc_number)
        calc_number += 1
        return calc

    _ = run_perf(geom, calc_getter, pals=1, mems=500, repeat=1, kind=kind)
