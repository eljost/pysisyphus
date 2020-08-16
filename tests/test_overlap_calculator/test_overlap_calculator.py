import pytest

import pysisyphus
from pysisyphus.calculators.OverlapCalculator import OverlapCalculator


@pytest.mark.parametrize(
    "cdds, fallback, no_mwfn, no_jmol", [
        (None    , None  , False, False),
        ("calc"  , None  , True , False),
        ("render", "calc", False, True),
        ("render", None  , True , True),
    ]
)
def test_cdds_fallback(cdds, fallback, no_mwfn, no_jmol, monkeypatch):
    # Disable Mwfn/Jmol as requested
    def mock_available(calculator):
        print("in mock!")
        if (calculator == "mwfn") and no_mwfn:
            print("mock nowfn")
            return False
        if (calculator == "jmol") and no_jmol:
            print("mock jmol")
            return False
        return True
    monkeypatch.setattr(pysisyphus.calculators.OverlapCalculator,
                        "available", mock_available)

    calc_kwargs = {
        "cdds": cdds,
    }
    calc = OverlapCalculator(**calc_kwargs)

    assert calc.cdds == fallback
