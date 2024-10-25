# [1] http://dx.doi.org/10.1016/j.cplett.2010.05.036
#     One-dimensional anharmonic oscillator: Quantum versus classical vibrational
#     partition functions
#     Beste, 2010

import numpy as np
import scipy as sp

from pysisyphus.constants import KB, KBAU, PLANCKAU, PLANCK


Polynomial = np.polynomial.Polynomial


_1OVERPLANCKAU = 1 / PLANCKAU
_2PI = 2 * np.pi


def classical_particle_density(
    x: float, poly: Polynomial, temperature: float, red_mass_au: float
) -> float:
    """Classical particle density.

    Eq. (9) in [1].

    Parameters
    ----------
    x
        Coordinate x at which the  provided polynomial is evaluated.
    poly
        1d-polynomial of arbitrary order used to evaluate the potential at
        the given coordinates in the numerical integration.
    temperature
        Temperature in Kelvin.
    red_mass_au
        Reduced mass of the oscillator in atomic units (NOT in amu!).


    Returns
    -------
    rho
        Classical particle density.
    """
    kbT = KBAU * temperature
    energy = poly(x) - poly(0.0)
    rho = _1OVERPLANCKAU * np.sqrt(_2PI * red_mass_au * kbT) * np.exp(-energy / kbT)
    return rho


def get_A1(x, deriv1: Polynomial, deriv2: Polynomial, temperature: float):
    """Correction term A1 in the Wigner-Kirkwood partition function.

    Eq. (6) in [1].

    Parameters
    ----------
    x
        Coordinate for which A2 is evaluated.
    deriv1
        First derivative of the potential w.r.t. coordinate x.
    deriv2
        Second derivative of the potential w.r.t. coordinate x.
    temperature
        Temperature in Kelvin.

    Returns
    -------
    A1
        Coefficient A1.
    """
    kbT = KBAU * temperature
    A1 = 1.0 / (12.0 * kbT) * deriv1(x) ** 2 - 1.0 / 6.0 * deriv2(x)
    return A1


def get_A2(
    x: float,
    deriv1: Polynomial,
    deriv2: Polynomial,
    deriv3: Polynomial,
    deriv4: Polynomial,
    temperature: float,
) -> float:
    """Correction term A2 in the Wigner-Kirkwood partition function.

    Eq. (7) in [1].

    Parameters
    ----------
    x
        Coordinate for which A2 is evaluated.
    deriv1
        First derivative of the potential w.r.t. coordinate x.
    deriv2
        Second derivative of the potential w.r.t. coordinate x.
    deriv3
        Third derivative of the potential w.r.t. coordinate x.
    deriv4
        Quartic derivative of the potential w.r.t. coordinate x.
    temperature
        Temperature in Kelvin.

    Returns
    -------
    A2
        Coefficient A2.
    """
    kbT = KBAU * temperature
    kbT2 = kbT**2
    deriv1_x = deriv1(x)
    deriv1_sq = deriv1_x**2
    deriv1_quart = deriv1_sq**2
    deriv2_x = deriv2(x)
    deriv2_sq = deriv2_x**2

    A2 = (
        (1.0 / (288.0 * kbT2) * deriv1_quart)
        - (11.0 / (360.0 * kbT) * deriv1_sq * deriv2_x)
        + (1.0 / 40.0 * deriv2_sq)
        + (1.0 / 30.0 * deriv1_x * deriv3(x))
        - (1.0 / 60.0 * kbT * deriv4(x))
    )
    return A2


def wigner_kirkwood_partfunc(
    poly: Polynomial,
    temperature: float,
    red_mass_au: float,
) -> tuple[float, float]:
    """Wigner-Kirkwood partition function (WKPF).

    Eqs. (5) to (7) and eq. (9) in [1].
    For small temperatures, e.g. 50 K and wavenumbers > 90 cm⁻¹ the WKPF
    is not reliable!

    Parameters
    ----------
    poly
        1d-polynomial of arbitrary order used to evaluate the potential at
        the given coordinates in the numerical integration.
    temperature
        Temperature in Kelvin.
    red_mass_au
        Reduced mass of the oscillator in atomic units (NOT in amu!).

    Returns
    -------
    q_wk
        Wigner-Kirkwood partition function.
    abserr
        Absolute error of the numerical integration.
    """
    kbT = KBAU * temperature
    lambda_factor = PLANCKAU**2 / (kbT**2 * 8.0 * red_mass_au * np.pi**2)
    lambda_factor2 = lambda_factor**2

    deriv1, deriv2, deriv3, deriv4 = [poly.deriv(m) for m in (1, 2, 3, 4)]

    def func(x):
        A1 = get_A1(x, deriv1, deriv2, temperature)
        A2 = get_A2(x, deriv1, deriv2, deriv3, deriv4, temperature)
        corr_factor = 1.0 + lambda_factor * A1 + lambda_factor2 * A2
        dens_classical = classical_particle_density(x, poly, temperature, red_mass_au)
        prod = dens_classical * corr_factor
        return prod

    y, abserr = sp.integrate.quad(func, -np.inf, np.inf)
    return y, abserr


def sos_partfunc_trunc(eigvals: np.ndarray, temperature: float, thresh=1e-4) -> float:
    """Truncated sum-over-states partition function.

    Eq. (2) in [1].

    Parameters
    ----------
    eigvals
        Energy-levels in Hartree.
    temperature
        Temperature in Kelvin.
    thresh
        Truncation threshold. When the contribution of an energy level falls below
        this threshold the calculation of the partition function is terminated.

    Returns
    -------
    q_sos
        Truncated sum-over-states partition function.
    """
    kbT = KBAU * temperature
    q_sos = 0.0
    for eigval in eigvals:
        dq = np.exp(-eigval / kbT)
        if abs(dq) <= thresh:
            break
        q_sos += dq
    return q_sos


def sos_partfunc(eigvals: np.ndarray, temperature: float) -> float:
    """Sum-over-states partition function.

    Eq. (2) in [1]. In contrast to [2] this partition function is not
    truncated but calculated from all provided eigenvalues. In the literature
    this approach is sometimes also referred to as "eigenvalue summation."

    Parameters
    ----------
    eigvals
        Energy-levels in Hartree.
    temperature
        Temperature in Kelvin.

    Returns
    -------
    q_sos
        Sum-over-states partition function.
    """
    q_sos = np.exp(-eigvals / (KBAU * temperature)).sum()
    return q_sos


def dln_sos_partfunc_dT(eigvals: np.ndarray, temperature: float) -> float:
    """Derivative of the sum-over-states partition function w.r.t. temperature.

    d ln(q_aq) / dT.

    Parameters
    ----------
    eigvals
        Energy-levels in Hartree.
    temperature
        Temperature in Kelvin.


    Returns
    -------
    dlnqdT
        Derivative of the natural logarithm of the sum-over-states partition function
        w.r.t. the temperature.
    """
    exp = np.exp(-eigvals / (KBAU * temperature))
    numerator = (eigvals / (KBAU * temperature**2) * exp).sum()
    dlnqdT = numerator / exp.sum()
    return dlnqdT


def anharmonic_classic_partfunc(
    poly: np.polynomial.Polynomial, temperature: float, red_mass_au: float
) -> float:
    """Classical anharmonic partition function.

    Eq. (9) in [1].

    Parameters
    ----------
    poly
        1d-polynomial of arbitrary order used to evaluate the potential at
        the given coordinates in the numerical integration.
    temperature
        Temperature in Kelvin.
    red_mass_au
        Reduced mass of the oscillator in atomic units (NOT in amu!).

    Returns
    -------
    q_ac
        Classical anharmonic partition function.
    """
    kbT = KBAU * temperature
    prefact = _1OVERPLANCKAU * np.sqrt(_2PI * red_mass_au * kbT)
    en_eq = poly(0.0)

    def func(x):
        en = poly(x) - en_eq
        return np.exp(-en / kbT)

    integral, _ = sp.integrate.quad(func, -np.inf, np.inf)
    q_ac = prefact * integral
    return q_ac


def harmonic_classic_partfunc(freq_si: float, temperature: float) -> float:
    """Classical partition function for a harmonic oscillator.

    Eq. (12) in [1].

    Parameters
    ----------
    freq_si
        Frequency of the oscillator in s⁻¹.
    temperature
        Temperature in Kelvin.

    Returns
    -------
    q_hc
        Classical harmonic partition function.
    """
    q_hc = KB * temperature / (PLANCK * freq_si)
    return q_hc


def harmonic_quantum_partfunc(freq_si: float, temperature: float) -> float:
    """Quantum harmonic partition function w/ bottom-of-well reference.

    Eq. (11) in [1].

    Parameters
    ----------
    freq_si
        Frequency of the oscillator in s⁻¹.
    temperature
        Temperature in Kelvin.

    Returns
    -------
    q_hq
        Quantum harmonic partition function.
    """
    quot = PLANCK * freq_si / (KB * temperature)
    q_hc = np.exp(-quot / 2.0) / (1.0 - np.exp(-quot))
    return q_hc
