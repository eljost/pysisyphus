import numpy as np
import pytest

from pysisyphus.calculators import Psi4
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.testing import using


@pytest.fixture
def azetidine():
    geom = geom_loader("lib:azetidine_hf_321g_opt.xyz", coord_type="redund")
    psi4_kwargs = {
        "method": "hf",
        "basis": "3-21g",
        "pal": 2,
    }
    calc = Psi4(**psi4_kwargs)
    geom.set_calculator(calc)
    return geom



@using("psi4")
def test_energy(azetidine):
    energy = azetidine.energy
    assert energy == pytest.approx(-171.1173045397321744)


@using("psi4")
def test_forces(azetidine):
    forces = azetidine.forces
    norm = np.linalg.norm(forces)
    assert norm == pytest.approx(2.567425805927507e-05, abs=1e-5)


@using("psi4")
def test_hessian(azetidine):
    H = azetidine.mw_hessian
    H = azetidine.eckart_projection(H)
    w, v = np.linalg.eigh(H)
    nus = eigval_to_wavenumber(w)
    ref_nus = np.array((3280.600842, 3322.467577, 3716.999024))
    np.testing.assert_allclose(nus[-3:], ref_nus, atol=1e-4)


@using("psi4")
@pytest.mark.parametrize(
    "pcm, ref_energy", [
        (  "cpcm", -171.1231505412408751),
        ("iefpcm", -171.1230983139995772),
    ]
)
def test_pcm_energy(azetidine, pcm, ref_energy):
    calc = azetidine.calculator
    calc.solvent = "water"
    calc.pcm = pcm
    energy = azetidine.energy
    assert energy == pytest.approx(ref_energy)


@using("psi4")
def test_mp2_energy():
    geom = geom_loader("lib:h2o.xyz")
    calc_kwargs = {
        "pal": 2,
        "method": "mp2",
        "basis": "6-31G*",
        "to_set": {
            "freeze_core": True,
        },
    }
    calc = Psi4(**calc_kwargs)
    geom.set_calculator(calc)
    energy = geom.energy
    assert energy == pytest.approx(-76.1959046366989412)
