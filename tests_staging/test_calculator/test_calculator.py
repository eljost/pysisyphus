#!/usr/bin/env python3

from pathlib import Path

import pytest

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.Turbomole import Turbomole


@pytest.mark.skip
def test_keep_ccre():
    path = Path("../../tests_staging/test_turbo_butadien_td_opt/butadiene")
    calc = Turbomole(path)
    ccre_pattern = "__ccre*"

    pat, multi, key = calc.prepare_pattern(ccre_pattern)
    assert pat == "ccre*"
    assert multi
    assert key == "ccres"


if __name__ == "__main__":
    test_keep_ccre()
