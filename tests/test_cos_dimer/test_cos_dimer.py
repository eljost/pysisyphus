import pytest

from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_hcn_neb_dimer_irc():
    run_dict = {
        "preopt": {
            "max_cycles": 3,
        },
        "interpol": {
            "type": "idpp",
            "between": 3,
        },
        "cos": {
            "type": "neb",
        },
        "opt": {
            "type": "qm",
            "align": True,
            "max_cycles": 10,
        },
        "tsopt": {
            "type": "dimer",
            "thresh": "gau_tight",
            "do_hess": True,
        },
        "irc": {
            "type": "eulerpc",
            "rms_grad_thresh": 1e-3,
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
        },
        "xyz": ["lib:hcn.xyz", "lib:hcn_iso_ts.xyz", "lib:nhc.xyz"]
    }
    results = run_from_dict(run_dict)

    assert results.ts_opt.is_converged
    assert results.ts_opt.cur_cycle == 7
    assert results.ts_geom.energy == pytest.approx(-92.2460427)

    irc = results.irc
    assert irc.forward_is_converged
    assert irc.forward_cycle == 29
    assert irc.backward_is_converged
    assert irc.backward_cycle == 36
