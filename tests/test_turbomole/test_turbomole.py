import numpy as np
import pytest

from pysisyphus.calculators import Turbomole
from pysisyphus.testing import using
from pysisyphus.helpers import geom_from_library, eigval_to_wavenumber


@pytest.mark.skip
@using("turbomole")
def test_turbomole_hessian():
    geom = geom_from_library("h2o_bp86_def2svp_opt.xyz")

    turbo_kwargs = {
        "control_path": "/scratch/programme/pysisyphus/tests/test_turbomole/ref_",
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    mw_H = geom.mw_hessian
    proj_H = geom.eckart_projection(mw_H)
    w, v = np.linalg.eigh(proj_H)
    nus = eigval_to_wavenumber(nus)
    print("nus / cm⁻¹:", nus)

    ref_nus = np.array((1607.81, 3684.62, 3783.64))
    np.testing.assert_allclose(nus[:3], ref_nus, atol=1e-2)


@using("turbomole")
@pytest.mark.parametrize(
        "control_path, ref_energy", [
        # Ground state
        ("/scratch/programme/pysisyphus/tests/test_turbomole/ref_", -76.36357867674),
        # Excited state
        ("/scratch/programme/pysisyphus/tests/test_turbomole/ref_exc", -76.0926146085),
        # ricc2
        ("/scratch/programme/pysisyphus/tests/test_turbomole/ref_ricc2", -75.8716368247),
        ],
)
def test_h2o_energy(control_path, ref_energy):
    geom = geom_from_library("h2o_bp86_def2svp_opt.xyz")
    turbo_kwargs = {
        "control_path": control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    energy = geom.energy

    assert energy == pytest.approx(ref_energy)


@using("turbomole")
@pytest.mark.parametrize(
        "control_path, ref_energy, ref_force_norm", [
        # Ground state
        ("/scratch/programme/pysisyphus/tests/test_turbomole/ref_",
            -76.36357867674, 1.30342385e-5),
        # Excited state gradient, TDDFT
        ("/scratch/programme/pysisyphus/tests/test_turbomole/ref_exc",
            -76.0926146085, 0.16006233),
        # Excited state gradient, ricc2
        ("/scratch/programme/pysisyphus/tests/test_turbomole/ref_ricc2",
            -75.8716368247, 0.15925937),
        ],
)
def test_h2o_forces(control_path, ref_energy, ref_force_norm):
    geom = geom_from_library("h2o_bp86_def2svp_opt.xyz")
    turbo_kwargs = {
        "control_path": control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces
    energy = geom.energy

    norm = np.linalg.norm(forces)

    assert norm == pytest.approx(ref_force_norm)
    assert energy == pytest.approx(ref_energy)
