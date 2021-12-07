import pytest

from pysisyphus.calculators import ORCA
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using
from pysisyphus.thermo import (
    print_thermoanalysis,
    get_thermoanalysis_from_hess_h5,
)
from pysisyphus.run import run_from_dict


@using("thermoanalysis")
def test_thermoanalysis(this_dir):
    hess_fn = this_dir / "h2o_hessian.h5"
    thermo = get_thermoanalysis_from_hess_h5(this_dir / "h2o_hessian.h5")

    assert thermo.M == pytest.approx(18.01528)
    assert thermo.dG == pytest.approx(0.002279584)


@pytest.fixture
def hcn_geom():
    """Optimized at HF/STO-3G"""
    geom = geom_loader("lib:hcn_sto3g_freq_ref.xyz")
    return geom


@using("pyscf")
@using("thermoanalysis")
def test_get_thermoanalysis(hcn_geom):
    hcn_geom.set_calculator(PySCF(basis="sto3g", verbose=0, pal=2))
    thermo = hcn_geom.get_thermoanalysis()
    print_thermoanalysis(thermo)

    assert thermo.dG == pytest.approx(-0.002, abs=1e-5)


@using("orca")
@using("thermoanalysis")
def test_hcn_thermo(hcn_geom):
    hcn_geom.set_calculator(ORCA(keywords="HF sto-3g"))
    thermo = hcn_geom.get_thermoanalysis()
    print_thermoanalysis(thermo, geom=hcn_geom)

    assert thermo.dG == pytest.approx(-0.002, abs=1e-5)


@using("pyscf")
@using("thermoanalysis")
def test_opt_h2o_do_hess():
    run_dict = {
        "geom": {
            "type": "redund",
            "fn": "lib:h2o.xyz",
        },
        "calc": {
            "type": "pyscf",
            "basis": "sto3g",
            "pal": 2,
            "verbose": 0,
        },
        "opt": {
            "thresh": "gau",
            "do_hess": True,
            "T": 398.15,
        },
    }
    run_result = run_from_dict(run_dict)


@using("thermoanalysis")
def test_print_thermo(this_dir):
    thermo, geom = get_thermoanalysis_from_hess_h5(
        this_dir / "hcn_orca_b973c_hessian.h5", return_geom=True
    )
    print_thermoanalysis(thermo, geom=geom)
