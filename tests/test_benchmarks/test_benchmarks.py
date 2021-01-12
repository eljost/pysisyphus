import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import ORCA


S22_SIZE = 22
BAKER_SIZE = 30
BAKER_TS_SIZE = 25


@pytest.mark.parametrize(
    "name, ref_size", [
    ("s22", S22_SIZE),
    ("baker", BAKER_SIZE),
    ("baker_ts", BAKER_TS_SIZE),
    ]
)
def test_benchmark_geoms(name, ref_size):
    bm = Benchmark(name)
    geoms = [geom for geom in bm.geoms]
    assert len(geoms) == ref_size


def test_benchmark_iter():
    bm = Benchmark("s22")
    assert len(list(bm)) == S22_SIZE


def test_calc_getter():
    def calc_getter(charge, mult):
        return ORCA(keywords="hf sto-3g", charge=charge, mult=mult)
    bm = Benchmark("baker_ts", calc_getter=calc_getter)
    geom = bm.get_geom(15)
    calc = geom.calculator
    assert calc is not None
    assert calc.charge == -1
    assert calc.mult == 1
    assert calc.keywords == "hf sto-3g"


def test_s22_exclude():
    exclude = (17, 18, 19)
    bm = Benchmark("s22", exclude=exclude)
    geoms = [geom for geom in bm.geoms]
    assert len(geoms) == S22_SIZE - len(exclude)
