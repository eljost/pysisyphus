import yaml
from pprint import pprint

import numpy as np
import pytest

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
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
            "climb_rms": 0.02,
        },
        "opt": {
            "type": "string",
            "stop_in_when_full": 3,
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
            # "basis": "321g",
            "basis": "sto3g",
            "verbose": 0,
        },
        "geom": {
            "type": "dlc",
            "fn": "lib:diels_alder_interpolated.trj",
        },
    }
    results = run_from_dict(run_dict)

    pprint(results)

    assert isinstance(results.cos, ChainOfStates)
    assert results.cos_opt.is_converged
    assert results.ts_opt.is_converged
    # assert results.ts_geom._energy == pytest.approx(-231.60321)  # 321g
    assert results.ts_geom._energy == pytest.approx(-230.036893)  # sto3g
    assert isinstance(results.ts_geom, Geometry)
    assert results.irc.forward_is_converged
    assert results.irc.backward_is_converged
    assert len(results.end_geoms) == 3
    assert results.calc_getter


@using("pyscf")
def test_run_results():
    run_dict = {
        "geom": {
            "type": "cart",
            "fn": "lib:h2o.xyz",
        },
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
                "verbose": 0,
            },
        },
        "geom": {
            "type": "cart",
            "fn": "lib:hcn_ts_hf_321g.xyz",
        },
    }
    results = run_from_dict(run_dict)

    assert results.opt.is_converged
    assert results.irc.backward_is_converged
    assert results.irc.forward_is_converged
    cnh, hcn = results.end_geoms
    cnh_ref_energy = -92.33966200907034
    hcn_ref_energy = -92.35408370090339
    # Sometimes mixed up, depending on the numpy version?!
    try:
        assert hcn.energy == pytest.approx(hcn_ref_energy)
        assert cnh.energy == pytest.approx(cnh_ref_energy)
    except AssertionError:
        assert hcn.energy == pytest.approx(cnh_ref_energy)
        assert cnh.energy == pytest.approx(hcn_ref_energy)


@using("pyscf")
def test_run_irc_constrained_endopt(this_dir):
    geom_path = this_dir / "ts_inp_run_constrained_endopt.xyz"
    ref_geom = geom_loader(geom_path)
    ref_c3d = ref_geom.coords3d.copy()

    constrain_ind = 0
    max_cycles = 5
    run_dict = {
        "geom": {
            "type": "cart",
            "fn": str(this_dir / "ts_inp_run_constrained_endopt.xyz"),
            "isotopes": [[constrain_ind, 1e9]],
        },
        "calc": {
            "type": "pyscf",
            "pal": 1,
            "basis": "321g",
            "verbose": 0,
        },
        "irc": {
            "type": "eulerpc",
            "max_cycles": max_cycles,
        },
        "endopt": {
            "geom": {
                "type": "redund",
                "coord_kwargs": {
                    "constrain_prims": [["atom", constrain_ind]],
                },
            },
            "max_cycles": max_cycles,
        },
    }
    results = run_from_dict(run_dict)
    feg, beg = results.end_geoms
    for end_geom in results.end_geoms:
        c3d = end_geom.coords3d
        np.testing.assert_allclose(c3d[constrain_ind], ref_c3d[constrain_ind])


@using("pyscf")
def test_new_style_yaml():
    run_dict = yaml.safe_load(
        """
    geom:
     type:
      tric:
       coord_kwargs:
        define_prims: [[BEND, 2, 0, 1]]
     fn: lib:h2o.xyz
    calc:
     type:
      pyscf:
       basis: sto3g
     pal: 6
    opt:
     type:
      rfo:
       hessian_init: fischer
       hessian_update: bfgs
     thresh: gau
     """
    )
    results = run_from_dict(run_dict)
    assert results.opt_geom.energy == pytest.approx(-74.965901183)
