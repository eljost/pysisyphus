import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.PreconSteepestDescent import PreconSteepestDescent
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators import Gaussian16, XTB
from pysisyphus.testing import using


@using("pyscf")
@pytest.mark.parametrize(
    "opt_cls, precon, ref_cycles",
    [
        (PreconSteepestDescent, True, 14),
        (PreconSteepestDescent, False, 54),
        (PreconLBFGS, True, 8),
        (PreconLBFGS, False, 7),
    ],
)
def test_water_hf_precon_opt(opt_cls, precon, ref_cycles):
    geom = geom_loader("lib:h2o_shaken.xyz")
    calc = PySCF(basis="sto3g")
    geom.set_calculator(calc)

    opt = opt_cls(geom, thresh="gau_tight", max_cycles=100, precon=precon)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycles
    assert geom.energy == pytest.approx(-74.96590119)


@using("gaussian16")
@pytest.mark.parametrize(
    "opt_cls, precon, ref_cycles",
    [
        (PreconSteepestDescent, True, 12),
        # (PreconSteepestDescent, False, 15),  # probably takes forever ...
        (PreconLBFGS, True, 6),
        (PreconLBFGS, False, 15),
    ],
)
def test_menthone(opt_cls, precon, ref_cycles):
    geom = geom_loader("lib:baker/29_menthone.xyz")
    calc = Gaussian16("PM6", pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "max_cycles": 100,
        "precon": precon,
        "overachieve_factor": 3.0,
    }

    opt = opt_cls(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycles


@pytest.mark.skip_ci
@using("xtb")
@pytest.mark.parametrize(
    "precon, precon_kind, ref_cycle",
    [
        (True, "full", 83),
        (True, "full_fast", 83),
        (True, "bonds", 93),
        (True, "bonds_bends", 83),
        # (False, None, 26),
    ],
)
def test_biaryl_precon(precon, precon_kind, ref_cycle):
    geom = geom_loader("lib:split.image_021.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "max_cycles": 100,
        "precon": precon,
        "precon_kind": precon_kind,
        "overachieve_factor": 3.0,
    }

    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    # Allow higher tolerance without preconditioner
    abs_ = 1e-4 if precon else 2e-3
    assert geom.energy == pytest.approx(-48.73588757, abs=abs_)


@using("pyscf")
@pytest.mark.parametrize(
    "precon_update, ref_cycles",
    [
        (1, 14),
        (5, 15),
        (10, 19),
    ],
)
def test_precon_update(precon_update, ref_cycles):
    geom = geom_loader("lib:h2o_shaken.xyz")
    calc = PySCF(basis="sto3g")
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau_tight",
        "max_cycles": 100,
        "precon_update": precon_update,
    }
    opt = PreconSteepestDescent(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycles
    assert geom.energy == pytest.approx(-74.96590119)


@using("pyscf")
def test_precon_freeze_atoms():
    geom_kwargs = {
        "coord_type": "cartesian",
        "freeze_atoms": [0, ],
    }
    geom = geom_loader("lib:h2o_shaken.xyz", **geom_kwargs)
    calc = PySCF(basis="sto3g")
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau",
        "precon_kind": "full",
    }
    opt = PreconLBFGS(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-74.96590115)
