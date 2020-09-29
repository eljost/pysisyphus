import pytest

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("xtb")
def test_endopt_barriers():
    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "splined_hei.xyz",
        },
        "tsopt": {
            "type": "rsirfo",
        },
        "irc": {
            "type": "eulerpc",
            "rms_grad_thresh": 0.015,
        },
        "endopt": {
            "thresh": "gau_loose",
        },
        "calc": {
            "type": "xtb",
            "pal": 1,
        },
    }
    results = run_from_dict(run_dict)

    assert results.ts_opt.is_converged
    assert results.ts_opt.cur_cycle == 1
    assert results.ts_geom.energy == pytest.approx(-5.38737424)
