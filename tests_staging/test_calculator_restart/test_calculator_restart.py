import numpy as np
import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators import ORCA, Gaussian16
from pysisyphus.testing import using


init_logging()


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, chk_ext",
    [
        pytest.param(
            ORCA, {"keywords": "HF def2-SV(P)", }, ".gbw",
            marks=using("orca")),
        pytest.param(
            Gaussian16, {"route": "HF/def2SVPP", }, ".fchk",
            marks=using("gaussian16")),
    ]
)
def test_restart(calc_cls, calc_kwargs, chk_ext):
    geom = geom_from_library("benzene.xyz")
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces

    ref_energy = -230.516816570
    assert geom.energy == pytest.approx(ref_energy)

    restart_info = calc.get_restart_info()
    chkfile = restart_info["chkfile"]
    assert chkfile is not None
    assert chkfile.endswith(chk_ext)

    calc2 = calc_cls(**calc_kwargs)
    calc2.set_restart_info(restart_info)
    geom2 = geom.copy()
    geom2.set_calculator(calc2)
    forces2 = geom2.forces

    np.testing.assert_allclose(forces2, forces, atol=1e-5)
