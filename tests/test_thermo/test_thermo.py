import pytest

from pysisyphus.calculators import ORCA
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.io import geom_from_hessian
from pysisyphus.testing import using
from pysisyphus.thermo import (
    print_thermoanalysis,
    get_thermoanalysis_from_hess_h5,
)
from pysisyphus.run import run_from_dict


@pytest.mark.skip
@using("thermoanalysis")
def test_thermoanalysis(this_dir):
    # H2O HF/321-G/RIJCOSX ORCA5
    hess_fn = this_dir / "h2o_hessian.h5"
    thermo = get_thermoanalysis_from_hess_h5(hess_fn, point_group="c2v")

    assert thermo.M == pytest.approx(18.01528)
    assert thermo.dG == pytest.approx(0.00412717, abs=1e-6)


@pytest.fixture
def hcn_geom():
    """Optimized at HF/STO-3G"""
    geom = geom_loader("lib:hcn_sto3g_freq_ref.xyz")
    return geom


@pytest.mark.skip
@using("pyscf")
@using("thermoanalysis")
def test_get_thermoanalysis(hcn_geom):
    hcn_geom.set_calculator(PySCF(basis="sto3g", verbose=4, pal=2))
    thermo = hcn_geom.get_thermoanalysis()
    print_thermoanalysis(thermo)

    assert thermo.dG == pytest.approx(-0.00029409, abs=1e-6)


@pytest.mark.skip
@using("orca")
@using("thermoanalysis")
def test_hcn_thermo(hcn_geom):
    hcn_geom.set_calculator(ORCA(keywords="HF sto-3g"))
    thermo = hcn_geom.get_thermoanalysis()
    print_thermoanalysis(thermo, geom=hcn_geom)

    assert thermo.dG == pytest.approx(-0.00029614, abs=1e-5)


@using("pyscf")
@using("thermoanalysis")
def test_opt_h2o_do_hess():
    T = 398.15
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
            "T": T,
        },
    }
    run_result = run_from_dict(run_dict)
    thermo = run_result.opt_geom.get_thermoanalysis(T=T)
    assert thermo.dG == pytest.approx(-0.00405113)


@using("thermoanalysis")
def test_print_thermo(this_dir):
    thermo, geom = get_thermoanalysis_from_hess_h5(
        this_dir / "hcn_orca_b973c_hessian.h5", return_geom=True
    )
    print_thermoanalysis(thermo, geom=geom)


@using("thermoanalysis")
@pytest.mark.parametrize(
    "id_, dG_ref", (
        # Ref values from ORCA logfiles
        (24, 0.62709533),
        (63, 0.62781152),
        (84, 0.62876245),
    )
)
def test_irc_h5(this_dir, id_, dG_ref):
    h5 = this_dir / f"irc_000.0{id_}.orca.h5"
    geom = geom_from_hessian(h5)
    thermo = geom.get_thermoanalysis()
    print_thermoanalysis(thermo)
    assert thermo.dG == pytest.approx(dG_ref, abs=2.5e-3)
