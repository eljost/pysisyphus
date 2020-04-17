import os
from pathlib import Path

import numpy as np
import pytest

from pysisyphus.calculators import ORCA, XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.modefollow import davidson
from pysisyphus.modefollow.NormalMode import NormalMode
from pysisyphus.testing import using


THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))
# init_logging()


@pytest.mark.parametrize(
    "precon", [
        True,
        False,
    ]
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_nu_pre, ref_cyc_pre, ref_nu, ref_cyc", [
        pytest.param(
            # 1791.60 cm⁻¹ from ORCA itself
            ORCA, {"keywords": "BP86 def2-SVP"},
            1793.274093, 3, 1793.293512, 7, marks=using("orca")),
        pytest.param(
            # 1691.888 cm⁻¹ from XTB itself
            XTB, {},
            1691.052689, 3, 1691.005160, 7, marks=using("xtb")),
        pytest.param(
            PySCF, {"basis": "321g",},
            # 1889.207638, 16, 1889.030654, 16, marks=using("pyscf")),
            # 1889.2027727, 16, 1889.012645, 16, marks=using("pyscf")),
            1889.2027727, 16, 1888.9849030442317, 16, marks=using("pyscf")),
    ],
)
def test_davidson_acet(precon, calc_cls, calc_kwargs,
                       ref_nu, ref_cyc, ref_nu_pre, ref_cyc_pre):
    geom = geom_loader("lib:acet_tm.xyz")

    calc_kwargs.update({
        "pal": 2,
    })
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    # Initial guess
    l = np.zeros_like(geom.coords).reshape(-1, 3)
    l[0][2] = 0.8
    l[1][2] = -0.6
    q = NormalMode(l.flatten(), geom.masses_rep)
    
    if precon:
        # Hessian for preconditioner
        # precon_calc = XTB(pal=2)
        # precon_geom = geom.copy()
        # precon_geom.set_calculator(precon_calc)
        # hessian_precon = precon_geom.eckart_projection(precon_geom.mw_hessian)
        # np.savetxt("hessian_precon", hessian_precon)
        hessian_precon = np.loadtxt("hessian_precon")
    else:
        hessian_precon = None

    davidson_kwargs = {
        "hessian_precon": hessian_precon,
    }
    result = davidson(geom, q, **davidson_kwargs)

    nu = result.nus[result.mode_ind]
    print(nu)

    ref_cycles = ref_cyc_pre if precon else ref_cyc
    nu_ref = ref_nu_pre if precon else ref_nu

    assert nu == pytest.approx(nu_ref, abs=0.5)
    assert result.cur_cycle+1 == ref_cycles
