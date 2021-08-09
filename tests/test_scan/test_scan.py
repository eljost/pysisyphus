import pytest

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
@pytest.mark.parametrize(
    "start, end, step_size", (
        (None, 4.0, None),  # Use initial value
        (3.0, 4.0, None),  # Start from 3.0
        (None, None, 0.2) # 3 steps with 0.2
    )
)
def test_h2o2_relaxed_scan(start, end, step_size):
    steps = 3
    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "lib:h2o2_hf_321g_opt.xyz",
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
            "verbose": 0,
        },
        "scan": {
            "type": "BOND",
            "indices": [2, 3],
            "start": start,
            "end": end,
            "steps": steps,
            "step_size": step_size,
            "opt": {
                "thresh": "gau",
            },
        },
    }
    results = run_from_dict(run_dict)
    # Original geometry is also returned
    assert len(results.scan_geoms) == (steps + 1)


@using("pyscf")
def test_h2o2_relaxed_scan_symmetric():
    steps = 3
    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "lib:h2o2_hf_321g_opt.xyz",
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
            "verbose": 0,
        },
        "scan": {
            "type": "BOND",
            "indices": [2, 3],
            "steps": 3,
            "step_size": 0.2,
            "symmetric": True,
            "opt": {
                "thresh": "gau",
            },
        },
    }
    results = run_from_dict(run_dict)
    # Original geometry is also returned
    assert len(results.scan_geoms) == 7
