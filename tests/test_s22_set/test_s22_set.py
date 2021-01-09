import pytest

from pysisyphus.benchmark_sets import get_s22_fns
from pysisyphus.calculators import ORCA
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using

import shutil


@using("orca")
@pytest.mark.parametrize(
    "fn", get_s22_fns()
)
def test_s22_set(fn):
    print(f"@Running: {fn}")
    geom = geom_loader(f"lib:s22/{fn}", coord_type="redund")

    id_ = int(fn[:2])
    calc_kwargs = {
        "keywords": "RI-MP2 6-31G** def2-SVP/C tightscf",
        "pal": 8,
        "mem": 1500,
        "charge": 0,
        "mult": 1,
    }
    calc = ORCA(**calc_kwargs)
    geom.set_calculator(calc)

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
