import numpy as np
import pytest

from pysisyphus.calculators import ORCA5, ConicalIntersection
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


init_logging()


@pytest.fixture
def calc():
    # Ground state
    gs_kwargs = {
        "keywords": "b3lyp def2-svp",
        "pal": 6,
        "calc_number": 0,
    }
    calc1 = ORCA5(**gs_kwargs)
    # S_1 state
    es_kwargs = gs_kwargs.copy()
    es_kwargs.update(
        {
            "numfreq": True,
            "blocks": "%tddft iroot 1 nroots 1 end",
            "calc_number": 1,
        }
    )
    calc2 = ORCA5(**es_kwargs)
    calc = ConicalIntersection(calc1, calc2)
    return calc


@using("orca5")
@pytest.mark.parametrize(
    "opt_cls, opt_kwargs, ref_energy",
    [
        (ConjugateGradient, {}, -78.25693885),
        (
            RFOptimizer,
            {
                "thresh": "gau",
                "trust_update": False,
            },
            -78.2487951,
        ),
    ],
)
def test_ci_opt(opt_cls, opt_kwargs, ref_energy, calc, this_dir):
    geom = geom_loader(this_dir / "ethene_init.xyz", coord_type="redund")
    geom.set_calculator(calc)
    opt = opt_cls(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(ref_energy)
