import os
from pathlib import Path
import shutil
import sys

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import ORCA
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


S22Bm = Benchmark("s22", coord_type="redund")


@pytest.mark.benchmark
@using("orca")
@pytest.mark.parametrize("fn, geom, charge, mult, ref_energy", S22Bm.geom_iter)
def test_s22_set(fn, geom, charge, mult, ref_energy, this_dir):
    calc = ORCA(
        keywords="RI-MP2 6-31G** def2-SVP/C tightscf",
        pal=8,
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
    backup_dir  = this_dir / "backup_gdiis_no_direction_test"
    if not backup_dir.exists():
        os.mkdir(backup_dir)

    def keep(fn, rm=False):
        shutil.copy(fn, backup_dir / f"{fn}.{id_:02}")
        if rm:
            os.remove(fn)

    keep("optimizer.log")
    keep("optimization.trj")
    keep("final_geometry.xyz")
    keep("optimization.h5")

    assert opt.is_converged
