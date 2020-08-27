from pprint import pprint

import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("pyscf")
def test_diels_alder_growing_string():
    run_dict = {
        "preopt": {
            "max_cycles": 5,
        },
        "cos": {
            "type": "gs",
            "fix_ends": True,
            "max_nodes": 8,
            "reparam_check": "rms",
            "climb": True,
            "climb_rms": 0.008,
        },
        "opt": {
            "type": "string",
            "stop_in_when_full": 0,
        },
        "tsopt": {
            "type": "rsirfo",
            "do_hess": True,
            "hessian_recalc": 5,
        },
        "irc": {
            "type": "eulerpc",
            "corr_func": "scipy",
            "rms_grad_thresh": 2.5e-3,
        },
        "endopt": {
            "fragments": True,
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
        },
        "xyz": "lib:diels_alder_interpolated.trj",
        "coord_type": "dlc",
    }
    results = run_from_dict(run_dict)

    pprint(results)

    assert len(results.preopt_xyz) == 3
    assert isinstance(results.cos, ChainOfStates)
    assert results.cos_opt.is_converged
    assert results.ts_opt.is_converged
    assert results.ts_geom._energy == pytest.approx(-231.6032000849316)
    assert isinstance(results.ts_geom, Geometry)
    assert results.irc.forward_is_converged
    assert results.irc.backward_is_converged
    assert len(results.end_geoms) == 3
    assert results.calc_getter


@using("pyscf")
def test_run_results():
    run_dict = {
        "xyz": "lib:h2o.xyz",
        "calc": {"type": "pyscf", "basis": "sto3g"},
    }
    results = run_from_dict(run_dict)

    calced_geoms = results.calced_geoms
    assert len(calced_geoms) == 1
    assert results.calc_getter


@using("pyscf")
def test_run_dimer_irc():
    """Quick test to see if the Dimer method works well with
    the rest."""
    run_dict = {
        "opt": {
            "type": "plbfgs",
            "do_hess": True,
        },
        "irc": {
            "type": "imk",
            "rms_grad_thresh": 0.001,
        },
        "endopt": {
            "fragments": True,
        },
        "calc": {
            "type": "dimer",
            "calc": {
                "type": "pyscf",
                "basis": "321g",
                "pal": 2,
            }
        },
        "xyz": "lib:hcn_ts_hf_321g.xyz",
    }
    results = run_from_dict(run_dict)

    assert results.opt.is_converged
    assert results.irc.backward_is_converged
    assert results.irc.forward_is_converged
    hcn, cnh = results.end_geoms
    assert hcn.energy == pytest.approx(-92.33966200907034)
    assert cnh.energy == pytest.approx(-92.35408370090339)
