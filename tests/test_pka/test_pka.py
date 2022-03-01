import pytest

from pysisyphus.drivers.pka import direct_cycle
from pysisyphus.testing import using


@using("thermoanalysis")
def test_direct_cycle(this_dir):
    acid_h5 = this_dir / "formic_opt_hessian.h5"
    base_h5 = this_dir / "formic_base_opt_hessian.h5"
    acid_solv_en = -189.698911494139
    base_solv_en = -189.233581882742
    pka = direct_cycle(acid_h5, base_h5, acid_solv_en, base_solv_en)
    assert pka == pytest.approx(9.7573997, abs=1e-6)
