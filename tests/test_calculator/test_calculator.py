from pathlib import Path

import numpy as np
import pytest

from pysisyphus.helpers_pure import results_to_json, json_to_results
from pysisyphus.run import run_from_dict
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


@using("pyscf")
@pytest.mark.parametrize(
    "run_func", (
        None,
        "get_energy",
        "get_forces",
        "get_hessian",
    )
)
def test_calculator_dump(run_func):
    run_dict = {
        "geom": {
            "fn": "lib:h2o.xyz",
        },
        "calc": {
            "type": "pyscf",
            "basis": "sto3g",
            "run_func": run_func,
        },
    }
    _ = run_from_dict(run_dict)
    assert Path("calculator_000.000.results").exists()
