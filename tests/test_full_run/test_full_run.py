from pprint import pprint

import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_full_run_hf_abstraction():
    run_dict = {
        "preopt": {
            "max_cycles": 25,
            "trust_max": 0.3,
            # "geom": {
                # "type": "tric",
            # },
        },
        "interpol": {
            "type": "redund",
            "between": 10,
        },
        "cos": {
            "type": "neb",
            "fix_ends": True,
        },
        "opt": {
            "type": "qm",
            "align": True,
            "max_cycles": 20,
        },
        "tsopt": {
            "type": "rsirfo",
            "do_hess": True,
        },
        "irc": {
            "type": "eulerpc",
            "corr_func": "scipy",
            "rms_grad_thresh": 5e-3,
        },
        "endopt": {
            "fragments": True,
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
        },
        "geom": {
            "type": "cart",
            "fn": ["lib:test_full_run_split.geom_008.xyz",
                  "lib:test_full_run_split.geom_034.xyz"
            ],
        },
    }
    results = run_from_dict(run_dict)

    assert results.preopt_first_geom
    assert results.preopt_last_geom
    assert isinstance(results.cos, ChainOfStates)
    assert results.ts_opt.is_converged
    assert results.ts_geom._energy == pytest.approx(-176.98452438)
    assert isinstance(results.ts_geom, Geometry)
    assert results.irc.forward_is_converged
    assert results.irc.backward_is_converged
    assert len(results.end_geoms) == 3
    assert results.calc_getter
