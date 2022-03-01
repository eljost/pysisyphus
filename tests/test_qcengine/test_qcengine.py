import numpy as np
import pytest

try:
    from pysisyphus.calculators.QCEngine import QCEngine
except ImportError:
    print("QCEngine import failed. Did you install it?")
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using
from pysisyphus.calculators import Turbomole


@using("turbomole")
@using("qcengine")
def test_qcengine_turbomole():
    geom = geom_loader("lib:h2o_guess.xyz")

    qce_kwargs = {
        "program": "turbomole",
        "model": {
            "method": "hf",
            "basis": "def2-SVP",
        },
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(-75.95615655854)
    assert norm == pytest.approx(0.11354367)


@pytest.mark.skip
@using("turbomole")
@using("qcengine")
def test_turbomole_hessian_compare(this_dir):
    geom = geom_loader("lib:h2o_bp86_def2svp_opt.xyz")

    qce_kwargs = {
        "program": "turbomole",
        "model": {
            "method": "b-p",
            "basis": "def2-SVP",
        },
        "keywords": {
            "ri": True,
            "grid": "m5",
        },
    }
    qce_calc = QCEngine(**qce_kwargs)
    geom.set_calculator(qce_calc)
    H = geom.hessian

    ref_geom = geom.copy()
    control_path = this_dir / "h2o_bp86_def2svp_control"
    ref_calc = Turbomole(control_path)
    ref_geom.set_calculator(ref_calc)
    H_ref = ref_geom.hessian

    assert geom.energy == pytest.approx(ref_geom.energy)
    np.testing.assert_allclose(H, H_ref, rtol=2e-4)


@using("mopac")
@using("qcengine")
def test_qcengine_mopac():
    geom = geom_loader("lib:h2o_guess.xyz")

    qce_kwargs = {
        "program": "mopac",
        "model": {
            "method": "pm6",
            "basis": None,
        },
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(-0.083446211)
    assert norm == pytest.approx(0.044034367)


@using("openmm")
@using("qcengine")
def test_qcengine_openmm():
    """
    conda install -c omnia -c conda-forge openmm
    conda install -c omnia -c conda-forge openforcefield
    """
    geom = geom_loader("lib:h2o_guess.xyz")

    qce_kwargs = {
        "program": "openmm",
        "model": {
            "method": "openmm",
            "basis": "openff-1.0.0",
            "offxml": "openff-1.0.0.offxml",
            "connectivity": ((0, 1, 1), (0, 2, 1)),
        },
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(6.4193809337e20)
    assert norm == pytest.approx(1.4649609864e22)

    # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # opt = RFOptimizer(geom)
    # opt.run()


@pytest.mark.skip
@using("gamess")
@using("qcengine")
def test_qcengine_gamess():
    geom = geom_loader("lib:h2o.xyz")

    qce_kwargs = {
        "program": "gamess",
        "model": {
            "method": "hf",
            "basis": "accd",
        },
        "keywords": {
            "contrl__ispher": 1,
        },
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(-76.0408384927)
    assert norm == pytest.approx(0.03196142235784051)
