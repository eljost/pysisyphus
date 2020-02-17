import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.PreconSteepestDescent import PreconSteepestDescent
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators import Gaussian16
from pysisyphus.testing import using


@using("pyscf")
@pytest.mark.parametrize(
    "opt_cls, precon, ref_cycles",
    [
        (PreconSteepestDescent, True, 7),
        (PreconSteepestDescent, False, 15),
        (PreconLBFGS, True, 7),
        (PreconLBFGS, False, 7),
    ]
)
def test_water_hf_precon_opt(opt_cls, precon, ref_cycles):
    geom = geom_from_library("h2o_shaken.xyz")
    calc = PySCF(basis="sto3g")
    geom.set_calculator(calc)

    opt = opt_cls(geom, thresh="gau_tight", precon=precon)
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
    ]
)
def test_menthone(opt_cls, precon, ref_cycles):
    geom = geom_from_library("baker/menthone.xyz")
    calc = Gaussian16("PM6", pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "max_cycles": 100,
        "precon": precon,
        "overachieve_factor": 3.,
    }

    opt = opt_cls(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycles
