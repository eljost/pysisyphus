import os
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
    "run_func",
    (
        None,
        "get_energy",
        "get_forces",
        "get_hessian",
    ),
)
def test_calculator_dump(run_func):
    out_dir = Path(".")
    out_fn = "calculator_000.000.results"
    try:
        os.remove(out_dir / out_fn)
    except FileNotFoundError:
        pass

    run_dict = {
        "geom": {
            "fn": "lib:h2o.xyz",
        },
        "calc": {
            "type": "pyscf",
            "basis": "sto3g",
            "run_func": run_func,
            "out_dir": out_dir,
        },
    }
    _ = run_from_dict(run_dict)
    assert (out_dir / out_fn).exists()
