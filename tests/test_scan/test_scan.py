from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_h2o2_relaxed_scan():
    end = 4.0
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
            "start": 3.0,
            "end": end,
            "steps": 5,
            "opt": {
                "thresh": "gau",
            },
        },
    }
    results = run_from_dict(run_dict)
    assert len(results.scan_geoms) == 6
