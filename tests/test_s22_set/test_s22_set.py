import shutil

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import ORCA
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


def calc_getter(charge, mult):
    return ORCA(keywords="RI-MP2 6-31G** def2-SVP/C tightscf", pal=6, mem=1500,
            charge=charge, mult=mult)
S22Bm = Benchmark("s22", coord_type="redund", calc_getter=calc_getter)

@using("orca")
@pytest.mark.parametrize(
    "fn, geom, ref_energy", S22Bm
)
def test_s22_set(fn, geom, ref_energy, results_bag):
    id_ = int(fn[:2])
    opt_kwargs = {
        "thresh": "gau",
        "dump": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged

    def keep(fn):
        shutil.copy(fn, f"{fn}.{id_:02}")
    keep("optimizer.log")
    keep("optimization.trj")
    keep("optimization.h5")
