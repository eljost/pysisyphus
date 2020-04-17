import numpy as np
import pytest

from pysisyphus.calculators import ORCA, XTB
from pysisyphus.helpers import geom_loader
from pysisyphus.modefollow import davidson
from pysisyphus.modefollow.NormalMode import NormalMode
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_nu", [
        pytest.param(
            # 1791.60 cm⁻¹ from ORCA itself
            ORCA,
            {"keywords": "BP86 def2-SVP"}, 1793.293512, marks=using("orca")),
        pytest.param(
            # 1691.888 cm⁻¹ from XTB itself
            XTB,
            {}, 1691.005160, marks=using("xtb")),
    ],
)
def test_davidson_acet(calc_cls, calc_kwargs, ref_nu):
    geom = geom_loader("lib:acet_tm.xyz")

    calc_kwargs.update({
        "pal": 2,
    })
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    geom.hessian

    # Initial guess
    l = np.zeros_like(geom.coords).reshape(-1, 3)
    l[0][2] = 0.8
    l[1][2] = -0.6
    q = NormalMode(l.flatten(), geom.masses_rep)
    
    nus, mode_ind = davidson(geom, q)

    assert nus[mode_ind]  == pytest.approx(ref_nu, abs=0.1)
