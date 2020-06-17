import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.dynamics.helpers import get_mb_velocities_for_geom
from pysisyphus.helpers import geom_loader, eigval_to_wavenumber
from pysisyphus.testing import using


@pytest.fixture
def geom():
    geom = geom_loader("lib:so3hcl_diss_ts_opt.xyz")
    calc = XTB(pal=2)
    geom.set_calculator(calc)
    return geom


@using("xtb")
def test_xtb_energy(geom):
    energy = geom.energy
    assert energy == pytest.approx(-20.46736065092)


@using("xtb")
def test_xtb_forces(geom):
    forces = geom.forces
    norm = np.linalg.norm(forces)
    assert norm == pytest.approx(0.000092380076)


@using("xtb")
def test_xtb_hessian(geom):
    Hmw = geom.mw_hessian
    Hmw = geom.eckart_projection(Hmw)

    w, v = np.linalg.eigh(Hmw)
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-571.77084564)
    assert nus[-1] == pytest.approx(1554.94377536)


@using("xtb")
def test_run_md(geom):
    T = 298.15
    seed = 20182503
    energy_ref = -20.46920283843
    velocities = get_mb_velocities_for_geom(geom, T=T, seed=seed)
    calc = geom.calculator
    geoms = calc.run_md(geom.atoms, geom.cart_coords, t=200, dt=0.1,
                        velocities=velocities)

    assert len(geoms) == 200
    last_geom = geoms[-1]
    energy = calc.get_energy(last_geom.atoms, last_geom.cart_coords)["energy"]
    assert energy == pytest.approx(energy_ref)
