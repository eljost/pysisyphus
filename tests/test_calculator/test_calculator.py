import os
from pathlib import Path

import numpy as np

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import results_to_json, json_to_results
from pysisyphus.testing import using


def test_results_to_json():
    size = 12
    results = {
        "energy": -1000.0,
        "forces": np.random.rand(size),
        "hessian": np.random.rand(size, size),
    }
    json_ = results_to_json(results)
    results_ = json_to_results(json_)
    for key, val in results_.items():
        ref_val = results[key]
        np.testing.assert_allclose(val, ref_val)


@using("xtb")
def test_calculator_dump():
    geom = geom_loader("lib:h2o.xyz")
    calc = XTB(dump=True)
    json_fn = calc.make_fn("results")
    try:
        os.remove(json_fn)
    except FileNotFoundError:
        pass
    geom.set_calculator(calc)
    geom.forces
    assert json_fn.exists()
