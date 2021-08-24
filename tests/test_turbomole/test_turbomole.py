# Tested against Turbomole 7.4.1

import numpy as np
from pathlib import Path
import pytest

from pysisyphus.calculators import Turbomole
from pysisyphus.testing import using
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber


@pytest.fixture
def geom():
    return geom_loader("lib:h2o_bp86_def2svp_opt.xyz")


@using("turbomole")
def test_turbomole_hessian(geom, this_dir):
    turbo_kwargs = {
        "control_path": this_dir / "./control_path_dft_gs",
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    mw_H = geom.mw_hessian
    proj_H = geom.eckart_projection(mw_H)
    w, v = np.linalg.eigh(proj_H)
    nus = eigval_to_wavenumber(w)
    print("nus / cm⁻¹:", nus)

    # # Turbomole reference values
    # turbo_ref_nus = np.array((1607.81, 3684.62, 3783.64))

    # I get slightly different values; probably from using different masses and
    # (slightly) different conversion factors when going from eigenvalues to
    # wavenumbers.
    ref_nus = np.array((1606.87, 3681.92, 3780.74))
    np.testing.assert_allclose(nus[-3:], ref_nus, atol=1e-2)


@using("turbomole")
@pytest.mark.parametrize(
        "control_path, ref_energy", [
        # Ground state
        ("./control_path_dft_gs", -76.36357867674),
        # Excited state
        ("./control_path_dft_es1", -76.0926146085),
        # ricc2
        ("./control_path_ricc2", -75.8716368247),
        ],
)
def test_h2o_energy(control_path, ref_energy, geom, this_dir):
    turbo_kwargs = {
        "control_path": this_dir / control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    energy = geom.energy

    assert energy == pytest.approx(ref_energy)


@using("turbomole")
@pytest.mark.parametrize(
        "control_path, ref_energy, ref_force_norm", [
        # Ground state
        ("./control_path_dft_gs",
            -76.36357867674, 1.30342385e-5),
        # Excited state gradient, TDDFT
        ("./control_path_dft_es1",
            -76.0926146085, 0.16006233),
        # Excited state gradient, ricc2
        ("./control_path_ricc2",
            -75.8716368247, 0.15925937),
        ],
)
def test_h2o_forces(control_path, ref_energy, ref_force_norm, geom, this_dir):
    turbo_kwargs = {
        "control_path": this_dir / control_path,
    }
    calc = Turbomole(**turbo_kwargs)
    geom.set_calculator(calc)

    forces = geom.forces
    energy = geom.energy

    norm = np.linalg.norm(forces)

    assert norm == pytest.approx(ref_force_norm, abs=1e-4)
    assert energy == pytest.approx(ref_energy)
