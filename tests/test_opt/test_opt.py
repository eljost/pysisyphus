#!/usr/bin/env python3

import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


@pytest.mark.parametrize(
    "calc, start, ref_cycle", [
    (AnaPot, (0.667, 1.609, 0.), 11)
    ]
)
def test_rfoptimizer(calc, start, ref_cycle):
    geom = calc.get_geom((0.667, 1.609, 0.))

    opt_kwargs = {
        "thresh": "gau_tight",
        "dump": False,
    }
    opt = RFOptimizer(geom)
    opt.run()
    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
