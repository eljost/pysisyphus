import numpy as np
import pytest

from pysisyphus.calculators import Gaussian16, Turbomole, ORCA
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_energy, ref_force_norm",
    [
        pytest.param(
            ORCA,
            {"keywords": "BP86 def2-SVP"},
            -40.473648820542,
            0.0539577447,
            marks=using("orca"),
        ),
        pytest.param(
            ORCA,
            # DoEq includes interaction between point charges and nuclear
            # charges of the actual atoms.
            {"keywords": "BP86 def2-SVP", "blocks": "%method doEQ true end"},
            -40.473658794066,
            0.0539577447,
            marks=using("orca"),
        ),
        pytest.param(
            Turbomole,
            {"control_path": "./methane_control_path"},
            -40.47560386,
            0.05402566536,
            marks=using("turbomole"),
        ),
        pytest.param(
            Gaussian16,
            {"route": "BP86 def2SVP"},
            -40.44843869845143,
            0.054086037,
            marks=using("gaussian16"),
        ),
        pytest.param(
            PySCF,
            # When check_mem is enabled the gibhut CI fails somehow...
            {"basis": "def2svp", "xc": "bp86", "check_mem": False},
            -40.473635092,
            0.05408810,
            marks=(using("pyscf"), pytest.mark.skip),
        ),
    ],
)
def test_with_point_charges(calc_cls, calc_kwargs, ref_energy, ref_force_norm, this_dir):
    geom = geom_loader("lib:methane_bp86_def2svp_opt.xyz")

    if "control_path" in calc_kwargs:
        calc_kwargs["control_path"] = this_dir / calc_kwargs["control_path"]

    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    point_charges = np.array(
        (
            (2.0, 1.78620256244410, 0.0, 0.3),
            (2.0, -1.78620256244410, 0.0, 0.3),
        )
    )

    results = calc.get_forces(geom.atoms, geom.coords, point_charges=point_charges)

    assert results["energy"] == pytest.approx(ref_energy)
    assert np.linalg.norm(results["forces"]) == pytest.approx(ref_force_norm, abs=1e-4)

    # results = calc.get_forces(geom.atoms, geom.coords, prepare_kwargs=None)
    # print(calc_cls, results["energy"], np.linalg.norm(results["forces"]))


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs, ref_energy, ref_charges",
    [
        pytest.param(
            Gaussian16,
            {"route": "BP86 def2SVP"},
            -40.4818752416,
            (-0.195353, 0.048838, 0.048838, 0.048838, 0.048838),
            marks=using("gaussian16"),
        ),
    ],
)
def test_parse_charges(calc_cls, calc_kwargs, ref_energy, ref_charges):
    geom = geom_loader("lib:methane_bp86_def2svp_opt.xyz")

    calc = calc_cls(**calc_kwargs)
    geom.set_calculator(calc)

    results = calc.get_forces(geom.atoms, geom.coords)

    assert results["energy"] == pytest.approx(ref_energy)

    charges = calc.parse_charges()
    np.testing.assert_allclose(charges, ref_charges, atol=1e-5)
