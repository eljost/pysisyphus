import itertools as it

import numpy as np
import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import XTB
from pysisyphus.irc.IRCDummy import IRCDummy
from pysisyphus.run import run_endopt
from pysisyphus.testing import using


@pytest.fixture
def ts():
    b = Benchmark("xtb_rx", calc_getter=calc_getter)
    geoms = b.get_geoms(0)
    ts_ = geoms[1]
    return ts_


def calc_getter(charge=0, mult=1):
    calc = XTB(pal=2, charge=charge, mult=mult, quiet=True)
    return calc


def gen_all_coords(ts):
    from pysisyphus.irc import EulerPC
    irc = EulerPC(ts, rms_grad_thresh=5e-4, hessian_recalc=10)
    irc.run()
    np.savetxt("all_coords", irc.all_coords)


@pytest.fixture
def irc(this_dir, ts):
    all_coords = np.loadtxt(this_dir / "all_coords")
    irc_ = IRCDummy(all_coords, ts.atoms)
    return irc_


@using("xtb")
@pytest.mark.parametrize(
    "fragments, geoms_expected", [
        (True, 2+1),
        (False, 1+1),
        ("total", 3+1),
    ]
)
def test_endopt(fragments, geoms_expected, ts, irc):
    geom = ts
    endopt_key = "rfo"
    endopt_kwargs = {
        "fragments": fragments,
        "geom": {
            "type": "redund",
        }
    }
    opt_results = run_endopt(irc, endopt_key, endopt_kwargs, calc_getter)

    flattened = list(it.chain(*opt_results))
    assert len(flattened) == geoms_expected
