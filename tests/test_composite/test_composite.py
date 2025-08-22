import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.Composite import Composite, get_GaMD_composite_calc
from pysisyphus.finite_diffs import finite_difference_gradient
from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


@using("orca")
def test_ch4_composite(this_dir):
    """G2MS/R from https://doi.org/10.1021/jp963019q"""

    geom = geom_loader(this_dir / "00_ch4.xyz")

    calc_kwargs = {
        "calcs": {
            "ccsdt": {
                "type": "orca",
                "keywords": "ccsd(t) 6-31G(d) tightscf",
            },
            "mp2_high": {
                "type": "orca",
                "keywords": "mp2 6-311+G(2df,2p) tightscf",
            },
            "mp2_low": {
                "type": "orca",
                "keywords": "mp2 6-31G(d) tightscf",
            },
        },
        "charge": 0,
        "mult": 1,
        "pal": 2,
        "final": "ccsdt + mp2_high - mp2_low",
    }
    calc = Composite(**calc_kwargs)
    geom.set_calculator(calc)
    en = geom.energy
    ref_energy = -40.42754204370099
    assert en == pytest.approx(ref_energy)


@using("pyscf")
def test_composite_run_dict(this_dir):
    run_dict = {
        "geom": {
            "type": "cart",
            "fn": str(this_dir / "00_ch4.xyz"),
        },
        "calc": {
            "type": "composite",
            "calcs": {
                "high": {
                    "type": "pyscf",
                    "basis": "321g",
                },
                "low": {
                    "type": "pyscf",
                    "basis": "sto3g",
                },
            },
            "final": "high - low",
            "pal": 2,
        },
    }
    results = run_from_dict(run_dict)
    geom = results.calced_geoms[0]
    assert geom._energy == pytest.approx(-0.250071439626311)


@using("xtb")
@pytest.mark.parametrize("dE", (-0.25, -0.1, 0.0, 0.1, 0.25, 0.5, 1.0))
def test_gamd_composite(dE):
    k = 0.12345
    E0 = -5.070431326355
    E = E0 + dE  # Boost threshold

    geom = geom_loader("lib:h2o.xyz")
    base = XTB(acc=0.0001)

    comp_calc = get_GaMD_composite_calc(base, k=k, E=E)
    geom.set_calculator(comp_calc)
    energy = geom.energy
    assert energy == pytest.approx(E0 + 0.5 * k * max(0.0, dE) ** 2)
    forces = geom.forces

    def scalar_func(coords):
        return geom.calculator.get_energy(geom.atoms, coords)["energy"]

    fd_grad = finite_difference_gradient(geom.cart_coords, scalar_func, step_size=1e-3)
    np.testing.assert_allclose(-fd_grad, forces, atol=2e-7)
