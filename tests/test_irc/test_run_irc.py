import pytest

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_from_library
from pysisyphus.run import run_irc
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "opt_ends, opt_geom_num", [
        pytest.param(False, 0, marks=using("xtb")),
        pytest.param(True, 2, marks=using("xtb")),
        pytest.param("fragments", 3, marks=using("xtb")),
    ]
)
def test_run_irc_opt_ends(opt_ends, opt_geom_num):
    geom = geom_from_library("hfabstraction_ts_opt_xtb.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    calc_ind = 0
    def calc_getter(index):
        nonlocal calc_ind
        calc = XTB(calc_number=calc_ind, pal=2)
        calc_ind += 1
        return calc

    irc_kwargs = {
        "type": "eulerpc",
        "opt_ends": opt_ends,
        "max_cycles": 15,

    }
    opt_geoms = run_irc(geom, irc_kwargs, calc_getter)
    assert len(opt_geoms) == opt_geom_num
