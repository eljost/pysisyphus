import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.testing import using


S22_SIZE = 22
BAKER_SIZE = 30
BAKER_TS_SIZE = 25
ZIMMERMAN_SIZE = 105
ZIMMERMAN_XTB_SIZE = 65
BIRKHOLZ_RX_SIZE = 20


@pytest.mark.parametrize(
    "name, ref_size",
    [
        ("s22", S22_SIZE),
        ("baker", BAKER_SIZE),
        ("baker_ts", BAKER_TS_SIZE),
        ("zimmerman", ZIMMERMAN_SIZE),
        ("zimmerman_xtb", ZIMMERMAN_XTB_SIZE),
        ("birkholz_rx", BIRKHOLZ_RX_SIZE),
    ],
)
def test_benchmark_geoms(name, ref_size):
    bm = Benchmark(name)
    assert len(list(bm)) == ref_size


@pytest.mark.parametrize("inv_exclude", [True, False])
def test_s22_exclude(inv_exclude):
    exclude = (17, 18, 19)
    bm = Benchmark("s22", exclude=exclude, inv_exclude=inv_exclude)
    ref_len = len(exclude) if inv_exclude else S22_SIZE - len(exclude)
    assert len(list(bm)) == ref_len


def test_s22_exclude():
    exclude = (17, 18, 19)
    bm = Benchmark("s22", exclude=exclude)
    assert len(list(bm)) == S22_SIZE - len(exclude)


@using("pyscf")
def test_calc_getter():
    def calc_getter(charge, mult):
        return PySCF(basis="sto3g", charge=charge, mult=mult)

    bm = Benchmark("baker_ts", calc_getter=calc_getter)
    geom = bm.get_geoms(15)
    calc = geom.calculator
    assert calc is not None
    assert calc.charge == -1
    assert calc.mult == 1
    assert calc.basis == "sto3g"
