import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.calculators import Composite
from pysisyphus.testing import using


@using("orca")
def test_ch4_composite(this_dir):
    """G2MS/R from https://doi.org/10.1021/jp963019q"""

    geom = geom_loader(this_dir / "00_ch4.xyz")

    calc_kwargs = {
        "from_dict": {
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
