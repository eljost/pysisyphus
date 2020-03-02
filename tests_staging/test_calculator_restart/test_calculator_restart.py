import numpy as np
import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators import ORCA
from pysisyphus.testing import using


init_logging()


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    [
        pytest.param(ORCA, {"keywords": "HF def2-SV(P)", }, marks=using("orca")),
    ]
)
def test_restart(calc_cls, calc_kwargs):
    geom = geom_from_library("benzene.xyz")
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces

    restart_info = calc.get_restart_info()
    assert restart_info["chkfile"] is not None

    calc2 = calc_cls(**calc_kwargs)
    calc2.set_restart_info(restart_info)
    geom2 = geom.copy()
    geom2.set_calculator(calc2)
    forces2 = geom2.forces

    np.testing.assert_allclose(forces2, forces, atol=1e-5)
