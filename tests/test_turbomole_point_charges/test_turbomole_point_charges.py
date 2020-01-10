import numpy as np
import pytest

from pysisyphus.calculators import Turbomole
from pysisyphus.helpers import geom_from_library
from pysisyphus.testing import using_turbomole


@using_turbomole
def test_turbomole_point_charges():
    geom = geom_from_library("methane_bp86_def2svp_opt.xyz")
    ctrl_path = "./methane_control_path"
    calc = Turbomole(ctrl_path)

    geom.set_calculator(calc)

    point_charges = np.array((
        (2.0,  1.78620256244410,   0.0,     0.3),
        (2.0, -1.78620256244410,   0.0,     0.3),
    ))

    prep_kwargs = {
        "point_charges": point_charges,
    }
    results = calc.get_forces(geom.atoms, geom.coords, prepare_kwargs=prep_kwargs)

    assert results["energy"] == pytest.approx(-40.47560386)
    assert np.linalg.norm(results["forces"]) == pytest.approx(0.05402566536)
