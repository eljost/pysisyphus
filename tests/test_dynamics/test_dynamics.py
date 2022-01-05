import numpy as np
import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.constants import VELO2E
from pysisyphus.dynamics.helpers import (
    kinetic_energy_from_velocities,
    kinetic_energy_for_temperature,
    temperature_for_kinetic_energy,
    remove_com_velocity,
    scale_velocities_to_temperatue,
    unscaled_velocity_distribution,
    get_mb_velocities_for_geom,
)
from pysisyphus.dynamics import md
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


def test_kinetic_energy_from_velocities():
    atoms = 5
    masses = np.ones(atoms)
    velocities = np.ones((atoms, 3))

    E_kin = kinetic_energy_from_velocities(masses, velocities)
    assert E_kin == pytest.approx(1.5 * atoms * VELO2E)


def test_kinetic_energy_calculations():
    atoms = 5
    T = 298.15

    E_kin = kinetic_energy_for_temperature(atoms, T)
    T_calc = temperature_for_kinetic_energy(atoms, E_kin)
    assert T == pytest.approx(T_calc)


def test_remove_com_velocity():
    atoms = 5
    velocities = np.ones((atoms, 3))
    masses = np.ones(atoms)
    v_ref = np.zeros_like(velocities)

    v = remove_com_velocity(velocities, masses, keep_norm=False)
    np.testing.assert_allclose(v, v_ref)


def test_scale_velocities_to_temperatue():
    atoms = 5
    masses = np.ones(atoms)
    T = 298.15
    E_ref = kinetic_energy_for_temperature(atoms, T)  # in Bohr
    v_unscaled = 2 * (np.random.rand(atoms, 3) - 0.5)
    v_scaled = scale_velocities_to_temperatue(masses, v_unscaled, T)
    E_scaled = kinetic_energy_from_velocities(masses, v_scaled)

    assert E_scaled == pytest.approx(E_ref)


def test_unscaled_velocity_distribution():
    atoms = 5
    masses = np.ones(atoms)
    T = 298.15
    seed = 20182503

    v = unscaled_velocity_distribution(masses, T, seed)
    v_ref = np.array(
        [
            [-17.5130082, 24.56713378, -3.46233695],
            [-1.12100338, -3.29778854, 15.19174986],
            [0.78607979, -31.36839825, -13.09479587],
            [23.69573447, 24.69377977, -8.79120253],
            [-16.84579461, 17.65313333, 13.93684248],
        ]
    )
    assert np.abs((v - v_ref)).sum() == pytest.approx(0.0, abs=1e-7)


@pytest.mark.parametrize(
    "remove_com_v, remove_rot_v",
    [
        (False, False),
        (True, True),
    ],
)
def test_get_mb_velocities_for_geom(remove_com_v, remove_rot_v):
    geom = geom_loader("lib:benzene.xyz")
    T = 298.15

    fixed_dof = 0
    if remove_com_v:
        fixed_dof += 3
    if remove_rot_v:
        fixed_dof += 3

    v = get_mb_velocities_for_geom(
        geom, T, remove_com_v=remove_com_v, remove_rot_v=remove_rot_v
    )  # in Bohr/fs
    E_kin = kinetic_energy_from_velocities(geom.masses, v)
    E_ref = kinetic_energy_for_temperature(
        len(geom.atoms), T, fixed_dof=fixed_dof
    )  # in Bohr
    assert E_kin == pytest.approx(E_ref)


def test_fixed_dof_kinetic_energy():
    atoms = 5
    T = 298.15

    E_kin_ref = kinetic_energy_for_temperature(atoms, T)

    fully_fixed = 3 * atoms
    E_kin_fully_fixed = kinetic_energy_for_temperature(atoms, T, fixed_dof=fully_fixed)
    assert E_kin_fully_fixed == 0.0

    translation_fixed = 3
    T_trans_fixed = temperature_for_kinetic_energy(
        atoms, E_kin_ref, fixed_dof=translation_fixed
    )
    assert T_trans_fixed > T


@using("pyscf")
def test_mb_velocities():
    geom = geom_loader("lib:h2o.xyz")
    geom.set_calculator(PySCF(basis="sto3g"))

    # Preoptimization
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()
    print()

    T = 298.15
    seed = 20182503
    v0 = get_mb_velocities_for_geom(geom, T, seed=seed).flatten()
    steps = 100
    dt = 0.5
    res = md(geom, v0, steps, dt)
    assert dt * steps / 1000 == pytest.approx(res.t_ps)

    # from pysisyphus.xyzloader import coords_to_trj
    # coords_to_trj("md.trj", geom.atoms, res.coords)


@using("xtb")
def test_thermostat():
    geom = geom_loader("lib:dynamics/10_water.xyz")
    from pysisyphus.calculators.XTB import XTB

    geom.set_calculator(XTB(pal=2))

    opt = RFOptimizer(geom, thresh="gau_loose", max_cycles=25)
    opt.run()

    T = 298.15
    seed = 20182503
    v0 = get_mb_velocities_for_geom(geom, T, seed=seed).flatten()
    # steps = 6000
    steps = 250
    dt = 0.5
    md_kwargs = {
        "v0": v0,
        "steps": steps,
        "dt": dt,
        "remove_com_v": True,
        "thermostat": "csvr",
        "timecon": 25,
    }
    res = md(geom, **md_kwargs)
    # assert dt * steps / 1000 == pytest.approx(res.t)

    import matplotlib.pyplot as plt

    E_tot = res.E_tot
    E_tot -= E_tot.mean()
    from pysisyphus.constants import AU2KJPERMOL

    E_tot *= AU2KJPERMOL
    T = res.T
    fig, (ax0, ax1) = plt.subplots(nrows=2)
    ax0.plot(E_tot)
    ax0.axhline(E_tot.mean())
    ax0.set_title(r"$\Delta$E$_{tot}$ / kJ mol⁻¹")

    ax1.plot(T)
    T_mean = T.mean()
    print("T_mean", T_mean, "K")
    ax1.axhline(T_mean)
    ax1.set_title(f"T / K, avg. = {T_mean:.2f} K")

    plt.tight_layout()
    plt.show()

    fig.savefig("md.pdf")
