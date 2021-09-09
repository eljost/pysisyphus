import numpy as np
import pytest

from pysisyphus.calculators import ORCA5, ConicalIntersection
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.testing import using


init_logging()


@pytest.fixture
def calc():
    gs_kwargs = {
        "keywords": "b3lyp def2-svp",
        "pal": 6,
    }
    es_kwargs = gs_kwargs.copy()
    es_kwargs["numfreq"] = True
    es_kwargs["blocks"] = "%tddft iroot 1 nroots 2 end"
    calc1 = ORCA5(**gs_kwargs)
    calc2 = ORCA5(**es_kwargs)
    calc = ConicalIntersection(calc1, calc2)
    return calc


@using("orca5")
def test_ci_opt(calc, this_dir):
    geom = geom_loader(this_dir / "ethene_init.xyz", coord_type="redund")

    geom.set_calculator(calc)

    opt = ConjugateGradient(geom, max_step=0.2, thresh="gau")
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-78.25705786)


@pytest.mark.skip
@using("orca5")
def test_ci_hessian(calc, this_dir):
    geom = geom_loader(this_dir / "ethene_b3lyp_def2svp_ci01_opt.xyz", coord_type="redund")
    geom.set_calculator(calc)

    H = geom.cart_hessian
    np.save("ethene_CI_hessian", H)
    w, v = np.linalg.eigh(H)
