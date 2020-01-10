import numpy as np
import pytest

from pysisyphus.calculators import Turbomole, ORCA
from pysisyphus.helpers import geom_from_library
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_energy, ref_force_norm",
    [
        pytest.param(
            ORCA,
            {"keywords": "BP86 def2-SVP"},
            -40.473648820542, 0.0539577447,
            marks=using("orca"),),
        # pytest.param(
            # ORCA,
            # {"keywords": "BP86 def2-SVP", "blocks": "%method doEQ true end"},
            # -40.448455709343, 0.0539577447,
            # marks=using("orca"),)
        pytest.param(
            Turbomole,
            {"control_path": "./methane_control_path"},
            -40.47560386, 0.05402566536,
            marks=using("turbomole"),)
])
def test_turbomole_point_charges(calc_cls, calc_kwargs, ref_energy, ref_force_norm):
    geom = geom_from_library("methane_bp86_def2svp_opt.xyz")

    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    point_charges = np.array((
        (2.0,  1.78620256244410,   0.0,     0.3),
        (2.0, -1.78620256244410,   0.0,     0.3),
    ))

    prep_kwargs = {
        "point_charges": point_charges,
    }
    results = calc.get_forces(geom.atoms, geom.coords, prepare_kwargs=prep_kwargs)

    assert results["energy"] == pytest.approx(ref_energy)
    assert np.linalg.norm(results["forces"]) == pytest.approx(ref_force_norm)
