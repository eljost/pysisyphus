import pytest

from pysisyphus.calculators import XTB
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import geom_loader
from pysisyphus.run import run_tsopt_from_cos
from pysisyphus.testing import using


@pytest.fixture
def cos(this_dir):
    geoms = geom_loader(this_dir / "input.cycle_019.trj")

    for i, geom in enumerate(geoms):
        calc = XTB(calc_number=i)
        geom.set_calculator(calc)
    cos = NEB(geoms)
    return cos


@using("xtb")
@pytest.mark.parametrize(
    "coord_type",
    [
        "cart",
        "dlc",
        "redund",
    ],
)
def test_run_tsopt_from_cos(coord_type, cos):
    calc_number = len(cos.images)

    def calc_getter():
        nonlocal calc_number

        calc = XTB(calc_number=calc_number)
        calc_number += 1
        return calc

    tsopt_key = "rsirfo"
    tsopt_kwargs = {
        "do_hess": True,
        "hessian_recalc": 5,
        "geom": {
            "type": coord_type,
            "coord_kwargs": {},
        },
    }
    opt_result = run_tsopt_from_cos(cos, tsopt_key, tsopt_kwargs, calc_getter)

    assert opt_result.geom.energy == pytest.approx(-11.44519302)


@using("xtb")
def test_run_tsopt_from_cos_dimer(cos):
    opt_result = run_tsopt_from_cos(
        cos,
        tsopt_key="dimer",
        tsopt_kwargs={
            "geom": {
                "type": "cart",
            }
        },
        calc_getter=XTB,
    )

    assert opt_result.geom.energy == pytest.approx(-11.44519302, abs=2e-5)
