import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.io.molden import get_xtb_nuc_charges
from pysisyphus.dynamics.helpers import get_mb_velocities_for_geom
from pysisyphus.elem_data import INV_ATOMIC_NUMBERS
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import eigval_to_wavenumber
from pysisyphus.testing import using
from pysisyphus.wavefunction.pop_analysis import mulliken_charges


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
    assert norm == pytest.approx(0.0000923641, abs=1e-6)


@using("xtb")
def test_xtb_hessian(geom):
    Hmw = geom.mw_hessian
    Hmw = geom.eckart_projection(Hmw)

    w, v = np.linalg.eigh(Hmw)
    nus = eigval_to_wavenumber(w)
    assert nus[0] == pytest.approx(-572.1710873)
    assert nus[-1] == pytest.approx(1556.0269008)


@using("xtb")
def test_run_md(geom):
    T = 298.15
    seed = 20182503
    energy_ref = -20.47232690901  # Obtained with xtb 6.4.0
    velocities = get_mb_velocities_for_geom(geom, T=T, seed=seed)
    calc = geom.calculator
    geoms = calc.run_md(
        geom.atoms, geom.cart_coords, t=200, dt=0.1, velocities=velocities
    )

    # trj = "\n".join([geom.as_xyz() for geom in geoms])
    # with open("md.trj", "w") as handle:
    # handle.write(trj)
    assert len(geoms) == 200
    last_geom = geoms[-1]
    energy = calc.get_energy(last_geom.atoms, last_geom.cart_coords)["energy"]
    assert energy == pytest.approx(energy_ref)


@using("xtb")
def test_parse_charges(geom):
    energy = geom.energy
    calc = geom.calculator
    charges = calc.parse_charges()
    print(charges)
    assert len(charges) == 6
    assert charges[0] == pytest.approx(1.55773)


def test_xtb_normal_termination(this_dir):
    assert XTB.check_termination(this_dir / "xtb_pass.out")
    assert not XTB.check_termination(this_dir / "xtb_crash.out")


@using("xtb")
def test_xtb_retry_calc(this_dir):
    geom = geom_loader(this_dir / "recover.xyz")
    geom.set_calculator(XTB(pal=6, retry_etemp=1000.0, retry_calc=1))
    en = geom.energy
    assert en == pytest.approx(-108.70462911348)


@using("xtb")
def test_xtb_stored_wavefunction():
    geom = geom_loader("lib:h2o.xyz")
    calc = XTB(pal=1, wavefunction_dump=True)
    geom.set_calculator(calc)
    geom.energy
    # Wavefunction already does some internal sanity checking
    wf = calc.get_stored_wavefunction()
    assert wf.charge == 0
    pa = mulliken_charges(wf)
    ref_charges = (-0.56619393, 0.28309696, 0.28309696)
    np.testing.assert_allclose(pa.charges, ref_charges, atol=1e-6)


def test_xtb_nuc_charges():
    # Parametrized up to Z = 86
    atomic_nums = np.arange(1, 87, dtype=int)
    all_atoms = [INV_ATOMIC_NUMBERS[Z] for Z in atomic_nums]
    mod_charges = get_xtb_nuc_charges(all_atoms)
    ecp_electrons = get_xtb_nuc_charges(all_atoms, as_ecp_electrons=True)
    np.testing.assert_allclose(mod_charges + ecp_electrons, atomic_nums)
