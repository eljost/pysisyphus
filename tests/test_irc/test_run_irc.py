from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_irc
from pysisyphus.testing import using


@using("xtb")
def test_run_irc():
    geom = geom_loader("lib:hfabstraction_ts_opt_xtb.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    def calc_getter(calc_number=0):
        calc = XTB(calc_number=calc_number, pal=2)
        return calc

    max_cycles = 20
    irc_kwargs = {
        "type": "eulerpc",
        "max_cycles": max_cycles,
    }
    irc_key = irc_kwargs.pop("type")
    irc = run_irc(geom, irc_key, irc_kwargs, calc_getter)
    # max_cycles in each direction + TS + 2 initial displacments
    ref_shape = (max_cycles + 1 + 1 + 1 + max_cycles, geom.coords.size)
    assert irc.all_coords.shape == ref_shape
