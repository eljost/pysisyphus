#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.tsoptimizers import *


@pytest.mark.parametrize(
    "opt_cls, ref_cur_cycle",
    [
        pytest.param(TRIM, 7),
        pytest.param(RSIRFOptimizer, 8),
        pytest.param(RSPRFOptimizer, 11),
])
def test_tshessian_opts(opt_cls, ref_cur_cycle):
    geom = AnaPot.get_geom((-0.6, 2.2, 0.))

    opt_kwargs = {
        "trust_radius": 0.2,
        "dump": False,
    }
    opt = opt_cls(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cur_cycle

    # calc = geom.calculator
    # calc.plot_opt(opt, show=True)
