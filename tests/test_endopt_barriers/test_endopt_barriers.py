import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using
from pysisyphus.run import run_from_dict


@pytest.fixture
def run_dict():
    return {
        "geom": {
            "fn": "lib:hcn_iso_ts_opt_xtb.xyz",
        },
        "irc": {
            "type": "eulerpc",
            "rms_grad_thresh": 0.015,
            "max_cycles": 5,
        },
        "endopt": {
            "thresh": "gau_loose",
            "max_cycles": 5,
        },
        "calc": {
            "type": "xtb",
            "pal": 1,
        },
    }


@using("thermoanalysis")
@using("xtb")
@pytest.mark.parametrize(
    "downhill, do_hess", [(True, True), (True, False), (False, True), (False, False)]
)
def test_hcn_endopt_barriers(run_dict, this_dir, downhill, do_hess):
    run_dict["endopt"]["do_hess"] = do_hess
    if downhill:
        run_dict["geom"]["fn"] = str(this_dir / "test_quick_downhill.xyz")
        run_dict["irc"]["downhill"] = True
    results = run_from_dict(run_dict)


# @pytest.mark.skip
@using("xtb")
def test_total_endopt_barriers(run_dict):
    run_dict = {
        "geom": {
            "fn": "lib:xtb_rx/00_c2no2.trj[1]",
        },
        "calc": {
            "type": "xtb",
            "pal": 2,
        },
        "irc": {
            "max_cycles": 10,
        },
        "endopt": {"fragments": "total", "geom": {"type": "tric"}},
    }
    results = run_from_dict(run_dict)
