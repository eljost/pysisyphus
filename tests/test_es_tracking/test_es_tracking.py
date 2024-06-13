import pytest

from pysisyphus.calculators import Gaussian16, ORCA, Turbomole
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "ovlp_type",
    (
        "tden",
        # "top",
    ),
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(
            ORCA,
            {
                "keywords": "uhf sto-3g tightscf",
                "blocks": "%tddft tda true nroots 10 end",
            },
            marks=using("orca"),
        ),
        pytest.param(
            Gaussian16,
            {
                "route": "uhf sto-3g tda=(nstates=10)",
            },
            marks=using("gaussian16"),
        ),
        pytest.param(
            Turbomole,
            {
                "simple_input": {
                    "basis": "sto-3g hondo",
                    "uhf": None,
                    "scfconv": 8,
                    "scfiterlimit": 200,
                    "scfinstab": "ucis",
                    "soes": "all 10",
                    "exopt": 2,
                    "denconv": "1d-7",
                },
            },
            marks=using("turbomole"),
        ),
    ),
)
def test_ucis_es_tracking(calc_cls, calc_kwargs, ovlp_type):
    geom = geom_loader("lib:h2o_hf_sto3g_orca_opt.xyz", coord_type="redund")

    calc_kwargs.update(
        {
            "root": 2,
            "track": True,
            "charge": 0,
            "mult": 1,
            "ovlp_type": ovlp_type,
        }
    )
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert geom.energy == pytest.approx(-74.642649, abs=5e-5)
    # Root flips from 2 to 3 in cycle 1
    assert geom.calculator.root == 3


@pytest.mark.parametrize(
    "ovlp_type",
    (
        "tden",
        "top",
    ),
)
@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    (
        pytest.param(
            ORCA,
            {
                "keywords": "rhf sto-3g tightscf",
                "blocks": "%tddft tda true nroots 8 end",
            },
            marks=using("orca"),
        ),
        pytest.param(
            Gaussian16,
            {
                "route": "rhf sto-3g tda=(nstates=8)",
            },
            marks=using("gaussian16"),
        ),
        pytest.param(
            Turbomole,
            {
                "simple_input": {
                    "basis": "sto-3g hondo",
                    "scfconv": 8,
                    "scfiterlimit": 200,
                    "scfinstab": "ciss",
                    "soes": "all 8",
                    "exopt": 5,
                    "denconv": "1d-7",
                },
            },
            marks=using("turbomole"),
        ),
    ),
)
def test_rcis_es_tracking(calc_cls, calc_kwargs, ovlp_type):
    geom = geom_loader("lib:h2o_hf_sto3g_orca_opt.xyz", coord_type="redund")

    calc_kwargs.update(
        {
            "root": 5,
            "track": True,
            "charge": 0,
            "mult": 1,
            "ovlp_type": ovlp_type,
        }
    )
    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    opt_kwargs = {
        "thresh": "gau",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert geom.energy == pytest.approx(-74.350827, abs=5e-5)
    assert geom.calculator.root == 5
