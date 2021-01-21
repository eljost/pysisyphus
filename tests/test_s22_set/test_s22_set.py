from pathlib import Path
import shutil

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import ORCA
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


S22Bm = Benchmark("s22", coord_type="redund")


@using("orca")
@pytest.mark.parametrize("fn, geom, charge, mult, ref_energy", S22Bm.geom_iter)
def test_s22_set(fn, geom, charge, mult, ref_energy, results_bag):
    calc = ORCA(
        keywords="RI-MP2 6-31G** def2-SVP/C tightscf",
        pal=6,
        mem=1500,
        charge=charge,
        mult=mult,
    )
    geom.set_calculator(calc)
    opt_kwargs = {
        "thresh": "gau",
        "dump": True,
        "gdiis": True,
        "line_search": False,
        "gdiis_test_direction": False,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    id_ = int(fn[:2])
    def keep(fn):
        shutil.copy(fn, Path("backup_gdiis_no_direction_test") / f"{fn}.{id_:02}")

    keep("optimizer.log")
    keep("optimization.trj")
    keep("optimization.h5")

    assert opt.is_converged
