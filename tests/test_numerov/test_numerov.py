# [1] https://doi.org/10.1039/C6CP06698D
#     Pushing the limit for the grid-based treatment of
#     Schrödinger's equation: a sparse Numerov approach for one,
#     two and three dimensional quantum problems
#     Kuenzer, Soraru, Hofer, 2016

import matplotlib.pyplot as plt
import numpy as np
import pytest
import scipy as sp

import pysisyphus.numerov as numerov
from pysisyphus.constants import (
    ANG2BOHR,
    KCALPERMOLPERANG2,
    AU2KCALPERMOL,
    AMU2AU,
    NU2AU,
)


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
    # As given for cc-pVQZ in Table 6 of [1]
    assert fundamental == pytest.approx(4162.1, abs=6e-2)
    _1stovertone = nus[2] - nus[0]
    assert _1stovertone == pytest.approx(8089.5, abs=7e-2)
    print(f"\n{fundamental=: >10.2f} cm⁻¹, {_1stovertone=: >10.2f} cm⁻¹")

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
