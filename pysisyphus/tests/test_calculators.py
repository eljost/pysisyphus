import numpy as np
import pytest

from pysisyphus.constants import ANG2BOHR
from pysisyphus.calculators import Gaussian09, Gaussian16, ORCA, PySCF, XTB, Psi4

try:
    from pysisyphus.calculators.PySCF import PySCF
except ModuleNotFoundError:
    # Dummy class so we can keep the PySCF reference in the test_h2o2 test
    class PySCF:
        pass


try:
    from pysisyphus.calculators.QCEngine import QCEngine
except ImportError:
    # Dummy class so we can keep the QCEngine reference in the test_h2o2 test
    class QCEngine:
        pass


from pysisyphus.helpers import Geometry
from pysisyphus.testing import using


@pytest.fixture
def h2o2():
    atoms = "H H O O".split()
    coords = (
        np.array(
            """
        -0.35738652   -2.42732105    0.01869782
        -0.55575507    0.31681832    1.27512527
        -0.71404094   -1.41851105    0.01869782
        -0.20069873   -0.69255477    1.27610279
    """.strip().split(),
            dtype=float,
        )
        * ANG2BOHR
    )
    return Geometry(atoms, coords)


@pytest.mark.parametrize(
    "calc_cls, calc_kwargs",
    [
        pytest.param(
            Gaussian09,
            {
                "route": "hf sto-3g",
            },
            marks=using("gaussian09"),
        ),
        pytest.param(
            Gaussian16,
            {
                "route": "hf sto-3g",
            },
            marks=using("gaussian16"),
        ),
        pytest.param(
            ORCA,
            {
                "keywords": "hf sto-3g",
            },
            marks=using("orca"),
        ),
        pytest.param(PySCF, {"basis": "sto3g"}, marks=using("pyscf")),
        pytest.param(
            Psi4,
            {
                "method": "hf",
                "basis": "sto-3g",
                "to_set": {
                    "scf_type": "direct",
                },
            },
            marks=using("psi4"),
        ),
        pytest.param(
            QCEngine,
            {
                "program": "gamess",
                "model": {"method": "hf", "basis": "sto"},
                "keywords": {
                    "basis__ngauss": 3,
                },
            },
            marks=using("gamess"),
            id="QCEngine_Gamess",
        ),
        pytest.param(
            QCEngine,
            {
                "program": "turbomole",
                "model": {
                    "method": "hf",
                    "basis": "sto-3g hondo",
                },
            },
            marks=using("turbomole"),
            id="QCEngine_Turbomole",
        ),
    ],
)
def test_h2o2_hf(calc_cls, calc_kwargs, h2o2):
    calc = calc_cls(**calc_kwargs)
    h2o2.set_calculator(calc)

    forces = h2o2.forces
    forces_norm = np.linalg.norm(forces)
    energy = h2o2.energy

    assert energy == pytest.approx(-148.722366, abs=1e-5)
    assert forces_norm == pytest.approx(0.15586495, abs=1e-5)


def test_h2o2_xtb(h2o2):
    """Tested with XTB 6.3.2"""
    calc = XTB(gfn=2)
    h2o2.set_calculator(calc)

    forces = h2o2.forces
    forces_norm = np.linalg.norm(forces)
    energy = h2o2.energy

    assert energy == pytest.approx(-9.02701527453)
    assert forces_norm == pytest.approx(0.130953039, abs=1e-5)
