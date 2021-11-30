import pytest

from pysisyphus.helpers_pure import standard_state_corr
from pysisyphus.io.hessian import geom_from_hessian
from pysisyphus.run import run_from_dict
from pysisyphus.drivers.barriers import do_endopt_ts_barriers
from pysisyphus.testing import using


@using("xtb")
@pytest.mark.parametrize(
    "do_hess",
    (
        False,
        True,
    ),
)
def test_endopt_total(do_hess):
    # ed, ts, prod = geom_loader("lib:xtb_rx/07_dacp_eth.trj")
    # print(geoms)
    run_dict = {
        "geom": {"type": "cart", "fn": "lib:xtb_rx/07_dacp_eth.trj[1]"},
        "calc": {
            "type": "xtb",
            "pal": 6,
            "quiet": True,
        },
        "irc": {
            "type": "eulerpc",
            "hessian_init": "hess_init_irc.h5",
        },
        "endopt": {
            "fragments": "total",
            "do_hess": do_hess,
        },
        "barriers": {
            "solv_calc": {
                "type": "xtb",
                "pal": 6,
                # "quiet": True,
                "gbsa": "water",
            },
            "do_standard_state_corr": True,
        },
    }
    results = run_from_dict(run_dict)
    print(results)


def test_standard_state_corr():
    T = 298.15
    p = 101325
    dG = standard_state_corr(T=T, p=p)
    assert dG == pytest.approx(0.003018805)


from pathlib import Path


@pytest.fixture
def ts_geom(this_dir):
    return geom_from_hessian(this_dir / "inp_ts_final_hessian.h5")


@pytest.fixture
def left_geoms(this_dir):
    return [
        geom_from_hessian(this_dir / "inp_forward_end_frag000_final_hessian.h5"),
    ]


@pytest.fixture
def right_geoms(this_dir):
    return [
        geom_from_hessian(this_dir / fn)
        for fn in (
            "inp_backward_end_frag000_final_hessian.h5",
            "inp_backward_end_frag001_final_hessian.h5",
        )
    ]


def solv_calc_getter(*args):
    class DummyCalc:
        def get_energy(self, atoms, coords):
            return {"energy": 0.0}

    return DummyCalc()


@pytest.mark.parametrize("do_standard_state_corr", (True, False))
@pytest.mark.parametrize("solv_calc_getter", (None, solv_calc_getter))
@pytest.mark.parametrize("do_thermo", (True, False))
def test_do_endopt_ts_barriers(
    do_standard_state_corr, solv_calc_getter, do_thermo, left_geoms, right_geoms, ts_geom
):

    do_endopt_ts_barriers(
        ts_geom,
        left_geoms,
        right_geoms,
        do_standard_state_corr=do_standard_state_corr,
        solv_calc_getter=solv_calc_getter,
        do_thermo=do_thermo,
    )
