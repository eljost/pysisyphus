# [1] https://doi.org/10.1039/C6CP06698D
#     Pushing the limit for the grid-based treatment of
#     Schrödinger's equation: a sparse Numerov approach for one,
#     two and three dimensional quantum problems
#     Kuenzer, Soraru, Hofer, 2016
# [2] https://doi.org/10.1016/j.cplett.2019.04.016
#     A periodic Numerov approach applied to the
#     torsional tunneling splitting in hydrogen peroxide,
#     aliphatic alcohols and phenol
#     Kuenzer, Hofer, 2019
# [3] https://doi.org/10.1119/1.1538575
#     Visualization and measurement of quantum rotational dynamics
#     Dimeo, 2003


import matplotlib.pyplot as plt
import numpy as np
import pytest
import scipy as sp

import pysisyphus.numerov as numerov
from pysisyphus.constants import (
    AU2EV,
    AU2J,
    AMU2AU,
    ANG2BOHR,
    AU2KJPERMOL,
    AU2KCALPERMOL,
    AU2NU,
    BOHR2M,
    HBAR,
    KCALPERMOLPERANG2,
    KG2AU,
    NU2AU,
)

from pysisyphus.elem_data import MASS_DICT


MASS_H = MASS_DICT["h"] * AMU2AU
MASS_D = MASS_DICT["d"] * AMU2AU


def plot_numerov(d, xs, energies, w, v, scale=1.0, show=True):
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.plot(xs, energies, c="k", ls="--", label="Pot")

    for i, psi in enumerate(v.T):
        wi = w[i]
        abs_psi2 = np.abs(psi) ** 2
        norm = sp.integrate.simpson(abs_psi2, x=xs, dx=d)
        ax.axhline(wi, xmin=-wi, xmax=wi, c="grey")
        ax.plot(xs, w[i] + scale * abs_psi2, label=f"State {i}")
        print(f"{i=: >3d}:, {w[i]=: 12.8f}, |ψ|²={norm: >12.8f}")
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(0, 100.0)
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("ΔE")
    if show:
        plt.show()
    return fig, ax


@pytest.mark.parametrize("num", range(251, 951, 100))
def test_harmonic_oscillator(num):
    nstates = 5
    # Parameters as given in the paper
    mass = 1.0 * AMU2AU
    k = 1098.018
    x0 = 0.0

    # Convert from kcal mol⁻¹ Å⁻² to atomic units
    k_au = k * KCALPERMOLPERANG2

    # Analytical solution
    omega = np.sqrt(k_au / mass)
    nu = np.arange(nstates)
    w_ref = omega * (nu + 0.5)

    def energy_getter(i, x):
        return 0.5 * k_au * (x - x0) ** 2

    xs = np.linspace(-2.5, 2.5, num=num) * ANG2BOHR
    d = xs[1] - xs[0]
    w, _ = numerov.run(xs, energy_getter, mass, nstates)
    max_err = np.abs(w - w_ref).max()
    print()
    with np.printoptions(suppress=True, precision=6):
        print(f"\t       d: {d: >12.8f} au")
        print("\t      Ref:", w_ref)
        print("\t     This:", w)
        print(f"\tmax(|Δ|): {max_err: >12.8e}")
    np.testing.assert_allclose(w, w_ref, atol=1e-6)
    # energies = np.array([energy_getter(None, x) for x in xs])
    # en0 = energy_getter(None, 0.0)
    # w_kcal = (w + en0) * AU2KCALPERMOL
    # plot_numerov(d, xs, energies * AU2KCALPERMOL, w_kcal, v, scale=2.0)


@pytest.mark.parametrize("num", range(151, 951, 100))
def test_morse_potential(num):
    nstates = 5
    mass = 1.0 * AMU2AU
    A = 174.147  # kcal/mol
    alpha = 1.776  # Angstrom
    x0 = 0.0

    # Convert from kcal mol⁻¹ to atomic units
    A_au = A / AU2KCALPERMOL
    alpha_au = alpha / ANG2BOHR

    # Analytical solution
    omega = alpha_au * np.sqrt(2 * A_au / mass)
    nu = np.arange(nstates)
    w_ref = omega * (nu + 0.5) - omega**2 / (4 * A_au) * (nu + 0.5) ** 2

    def energy_getter(i, x):
        return A_au * (1.0 - np.exp(-alpha_au * (x - x0))) ** 2

    xs = np.linspace(-1.5, 2.5, num=num) * ANG2BOHR
    d = xs[1] - xs[0]
    w, _ = numerov.run(xs, energy_getter, mass, nstates)
    max_err = np.abs(w - w_ref).max()
    print()
    with np.printoptions(suppress=True, precision=6):
        print(f"\t      d: {d: >12.8f} au")
        print("\t     Ref:", w_ref)
        print("\t    This:", w)
        print(f"\tmax(|Δ|): {max_err: >12.8e}")
    np.testing.assert_allclose(w, w_ref, atol=1e-6)


def test_hydrogen(this_dir):
    """Hydrogen at FCI/cc-pVQZ level of theory."""
    data = np.loadtxt(this_dir / "03_h2_scan_qz.dat", skiprows=1)
    xs, energies = data.T
    xs *= ANG2BOHR
    d = xs[1] - xs[0]
    mass_h = 1.007825032 * AMU2AU
    energies -= energies.min()
    # Reduced mass of H2
    red_mass = mass_h * mass_h / (mass_h + mass_h)

    def energy_getter(i, x):
        return energies[i]

    w, v = numerov.run(xs, energy_getter, red_mass, nstates=15)

    nus = w / NU2AU
    fundamental = nus[1] - nus[0]
    _1stovertone = nus[2] - nus[0]
    print(f"\n{fundamental=: >10.2f} cm⁻¹, {_1stovertone=: >10.2f} cm⁻¹")
    # As given for cc-pVQZ in Table 6 of [1]
    assert fundamental == pytest.approx(4162.1, abs=6e-2)
    assert _1stovertone == pytest.approx(8089.5, abs=7e-2)

    """
    energies_kcal = energies * AU2KCALPERMOL
    w_kcal = w * AU2KCALPERMOL
    xs_ang = xs / ANG2BOHR

    scale = 2.0
    # Compare figure to Fig. (3) in [1]
    _, ax = plt.subplots(figsize=(16, 8))
    ax.plot(xs_ang, energies_kcal, c="k", ls="--", label="Pot")

    for i, psi in enumerate(v.T):
        abs_psi2 = np.abs(psi) ** 2
        norm = sp.integrate.simpson(abs_psi2, x=xs, dx=d)
        ax.plot(xs_ang, w_kcal[i] + scale * abs_psi2, label=f"State {i}")
        ax.axhline(w_kcal[i], c="k", ls=":", alpha=0.5)
        print(f"{i=: >3d}:, {w_kcal[i]=: 14.8f}, |ψ|²={norm: >12.8f}")
    ax.set_xlim(xs_ang[0], xs_ang[-1])
    ax.set_ylim(-5.0, 125.0)
    # ax.legend()
    ax.set_xlabel("x / Å")
    ax.set_ylabel("ΔE / kcal mol⁻¹")
    # plt.show()
    """


def red_mass(m1, m2):
    return m1 * m2 / (m1 + m2)


@pytest.mark.parametrize(
    "mass_1, mass_2, ref",
    (
        # Reference values from Table 1 in [2].
        (MASS_H, MASS_H, 10.50),
        (MASS_D, MASS_H, 4.72),
        (MASS_D, MASS_D, 1.22),
    ),
)
def test_h2o2_ref(mass_1, mass_2, ref, this_dir):
    # Data from the SI of [2]. First column is in Angstrom, second column is in
    # kcal/mol. The data in the 1st column was (probably) obtained vis the following
    # procedure:
    #   1.) Relaxed scan around the O-O bond w/ appropriate spacing, yielding
    #       torsion angles and energies.
    #   2.) To convert the torsion angles into a Cartesian coordinate the Cartesian
    #       coordinate differences between all involved atoms (in this case all four
    #       atoms) between succesive steps in the scan were calculated and accumulated.
    #       This yields something like a total Cartesian displacement TOT_DISPL of
    #       all atoms along the scan.
    #   3.) Finally, the calculated data was probably splined using the Cartesian
    #       differences and re-evalualted on an equidistant grid in the interval
    #       [0, TOT_DISPL] w/ the appropriate number of points.
    #       The splining is probably necessary as the Cartesian differences are
    #       probably not constant.
    pot = np.loadtxt(this_dir / "h2o2_pot.dat")
    xs, energies = pot.T
    xs *= ANG2BOHR
    energies /= AU2KCALPERMOL

    energies -= energies.min()

    def energy_getter(i, x):
        return energies[i]

    mass = red_mass(mass_1, mass_2)

    w, v = numerov.run(
        xs,
        energy_getter,
        mass,
        nstates=13,
        accuracy=10,
        periodic=True,
    )
    ts = (w[1] - w[0]) * AU2NU
    print(f"Tunnel splitting: {ts: >10.2f} cm⁻¹")
    assert ts == pytest.approx(ref, abs=6e-3)

    # fig, ax = plt.subplots()
    # energies_kcal = energies * AU2KCALPERMOL
    # w_kcal = w * AU2KCALPERMOL
    # ax.plot(xs, energies_kcal)

    # for i, state in enumerate(v.T):
    # color = "red" if i % 2 == 0 else "blue"
    # ax.axhline(w_kcal[i], c="k", ls=":")
    # ax.plot(xs, w_kcal[i] + (0.375 * state), color=color)

    # ax.set_xlabel("Torsion / rad")
    # ax.set_ylabel("ΔE / kcal mol⁻¹")
    # ax.set_title(f"tunnel splitting: {ts: >10.4} cm⁻¹")
    # plt.show()


CH3_INERTIA_SI = 5.3e-47  # kg m²
CH3_INERTIA_AU = CH3_INERTIA_SI / BOHR2M**2 * KG2AU
CH3_B_AU = HBAR**2 / 2.0 / CH3_INERTIA_SI / AU2J  # in au
# One full rotation/360°
# Drop last point of grid, to make it suitable for a periodic calculation
CH3_GRID = np.linspace(0, 2.0 * np.pi, 2 * 360 + 1)[:-1]


def test_ch3i_free_rotor():
    """Free CH3 rotor in CH3I. Example from [3]."""
    nstates = 9

    # Potential is zero
    def energy_getter(i, x):
        return 0.0

    # Convert from kg m² to moment of inertia in atomic units (mass_au * Bohr**2)
    w, v = numerov.run(
        CH3_GRID, energy_getter, CH3_INERTIA_AU, periodic=True, nstates=nstates
    )

    # Set up reference energy leves; all states for j > 0 are doubly degenerate
    # Rotational quantum number js
    js = np.repeat(np.arange((nstates - 1) // 2 + 1), 2)[1:]  # Drop first 0
    w_ref = CH3_B_AU * js**2
    np.testing.assert_allclose(w, w_ref, atol=1e-12)

    """
    scale = 3e-5
    fig, ax = plt.subplots()
    ax.axhline(0.0, label="Zero")
    for i, wi in enumerate(w):
        ax.axhline(wi, c="k", ls=":", label=f"state {i}")
        ax.axhline(w_ref[i], c="red", ls="-.", label=f"ref {i}")
        ax.plot(CH3_GRID, wi + scale * v[:, i])
    ax.set_xlim(0, CH3_GRID[-1])
    ax.set_ylim(-2e-5, 1.25 * w[-1])
    ax.set_xlabel("deg / rad")
    ax.set_ylabel("E / au")
    fig.tight_layout()
    plt.show()
    """


@pytest.mark.parametrize(
    "V3_eV, atol0, atol1",
    (
        (41e-3, 1.5e-5, 8.0e-5),
        (22.5e-3, 1.6e-5, 9.0e-5),
    ),
)
def test_ch3i_hindered_rotor(V3_eV: float, atol0: float, atol1: float):
    """Hindered CH3 rotor in solid CH3I. Example from [3]."""
    nstates = 3 * 3
    # V3 = 41.0  # mEV
    V3 = V3_eV / AU2EV
    # Eq. (5) in [3] for high barries
    omega_0 = 3 * np.sqrt(V3 / (2.0 * CH3_INERTIA_AU))

    # Potential term from eq. (3) in [3]
    energies = V3 / 2.0 * (1.0 - np.cos(3 * CH3_GRID))
    energies -= energies.min()

    def energy_getter(i, x):
        return energies[i]

    # Convert from kg m² to moment of inertia in atomic units (mass_au * Bohr**2)
    w, v = numerov.run(
        CH3_GRID, energy_getter, CH3_INERTIA_AU, periodic=True, nstates=nstates
    )
    print()
    for i, wi in enumerate(w):
        print(f"{i:03d}: w={wi*AU2KJPERMOL: >8.3} kJ mol⁻¹")

    nlib = np.repeat((0, 1), 3)
    w_ref = (0.5 + nlib) * omega_0
    # Here, the allowed absolute error is quite high, because the formulate
    # for the harmonic oscillator does not quite fit.
    # The first three states are correct to 1.5e-5 au
    np.testing.assert_allclose(w[:3], w_ref[:3], atol=atol0)
    np.testing.assert_allclose(w[3:6], w_ref[3:], atol=atol1)

    """
    scale = 1e-4
    fig, ax = plt.subplots()
    ax.plot(CH3_GRID, energies)
    for i, wi in enumerate(w):
        ax.axhline(wi, c="k", ls=":", label=f"state {i}")
        # ax.axhline(w_ref[i], c="red", ls="-.", label=f"ref {i}")
        ax.plot(CH3_GRID, wi + scale * v[:, i])
    # ax.set_ylim(-2e-5, 1.25 * w[-1])
    ax.set_xlim(0, CH3_GRID[-1])
    ax.set_ylim(0.0, 1.25 * w[-1])
    ax.set_xlabel("deg / rad")
    ax.set_ylabel("E / au")
    fig.tight_layout()
    plt.show()
    """
