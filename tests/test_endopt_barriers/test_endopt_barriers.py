import pytest

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
    "downhill, do_thermo", [(True, True), (True, False), (False, True), (False, False)]
)
def test_hcn_endopt_barriers(run_dict, this_dir, downhill, do_thermo):
    run_dict["endopt"]["do_thermo"] = do_thermo
    if downhill:
        run_dict["geom"]["fn"] = str(this_dir / "test_quick_downhill.xyz")
        run_dict["irc"]["downhill"] = True
    results = run_from_dict(run_dict)
