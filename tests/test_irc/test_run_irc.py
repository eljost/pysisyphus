import pytest

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_irc, run_endopt
from pysisyphus.testing import using


@using("xtb")
@pytest.mark.parametrize(
    "fragments, opt_geom_num", [
        (False, 2),
        (True, 3),
    ]
)
def test_run_irc_opt_ends(fragments, opt_geom_num):
    geom = geom_loader("lib:hfabstraction_ts_opt_xtb.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    def calc_getter(calc_number=0):
        calc = XTB(calc_number=calc_number, pal=2)
        return calc

    irc_kwargs = {
        "type": "eulerpc",
        "max_cycles": 20,

    }
    irc_key = irc_kwargs.pop("type")
    irc = run_irc(geom, irc_key, irc_kwargs, calc_getter)

    endopt_key = "rfo"
    endopt_kwargs = {
        "thresh": "gau_tight",
        "fragments": fragments,
        "do_hess": True,
        "geom": {
            "type": "redund",
        },
    }
    fw_results, bw_results, _ = run_endopt(irc, endopt_key, endopt_kwargs, calc_getter)
    import pdb; pdb.set_trace()  # fmt: skip
    assert len(fw_results) + len(bw_results) == opt_geom_num
