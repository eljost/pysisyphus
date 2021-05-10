import pytest

from pysisyphus.calculators import Composite
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
        "pal": 6,
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
                }
            },
            "final": "high - low",
            "pal": 4,
        }
    }
    results = run_from_dict(run_dict)
    geom = results.calced_geoms[0]
    assert geom._energy == pytest.approx(-0.250071439626311)
