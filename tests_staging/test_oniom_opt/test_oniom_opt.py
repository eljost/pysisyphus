import pytest

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("xtb")
def test_oniom_opt():
    run_dict = {
        "geom": {
            "type": "redund",
            # frag.py raffinose.xyz 21 14 33 10
            "fn": "lib:birkholz/raffinose.xyz",
        },
        "calc": {
            "type": "oniom",
            "calcs": {
                "real": {
                    "type": "xtb",
                    "gfn": 0,
                    "pal": 6,
                },
                "high": {
                    "type": "xtb",
                    "gfn": 2,
                    "pal": 6,
                },
            },
            "models": {
                "high": {
                    "inds": "11..21,33..36,50..57",
                    "calc": "high",
                },
            },
        },
        "opt": {
            "align": True,
            "thresh": "gau",
        },
    }
    results = run_from_dict(run_dict)

    assert results.opt_geom.energy == pytest.approx(-110.75489874)
