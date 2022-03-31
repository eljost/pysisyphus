# [1] https://doi.org/10.1002/ange.201914943
#     Heavy-Atom Tunneling in Organic Reactions
#     Castro, Karney, 2020
# [2] https://doi.org/10.1002/wcms.1165
#     Theory and simulation of atom tunneling in chemical reactions
#     Kästner, 2014
# [3] https://doi.org/10.1002/qua.25686
#     Eyringpy
#     Dzib, Merino et al, 2018
# [4] https://pubs.acs.org/doi/pdf/10.1021/j100809a040
#     TUNNELLING CORRECTIONS FOR UNSYMMETRICAL ECKART POTENTIAL ENERGY BARRIERS
#     Johnston, Heicklen, 1961
# [5] https://dx.doi.org/10.6028%2Fjres.086.014
#     A Method of Calculating Tunneling Corrections For Eckart Potential Barriers
#     Brown, 1981

import numpy as np
import numpy.typing as npt
import scipy.integrate as integrate

from pysisyphus.constants import KB, KBAU, PLANCK, PLANCKAU, AU2SEC


def eyring_rate(
    barrier_height: float,
    temperature: npt.ArrayLike,
    trans_coeff: float = 1.0,
) -> npt.NDArray[np.float64]:
    """Rate constant in 1/s from the Eyring equation.

    See https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456
    "Reaction Rates and Barriers" on p. 3234 and eq. (8).

    Parameters
    ----------
    barrier_height
        Barrier height (energy, enthalpy, gibbs energy, ...) in Hartree.
    temperature
        Temperature in Kelvin.
    trans_coeff
        Unitless transmission coefficient, e.g., obtained from Wigner or Eckart
        correction.

    Returns
    -------
    rate
        Eyring reaction rate in 1/s.
    """
    temperature = np.array(temperature, dtype=float)
    prefactor = trans_coeff * KB * temperature / PLANCK
    rate = prefactor * np.exp(-barrier_height / (KBAU * temperature))
    return rate


def harmonic_tst_rate(
    barrier_height: float,
    temperature: float,
    rs_part_func: float,
    ts_part_func: float,
    trans_coeff: float = 1.0,
) -> float:
    """Rate constant in 1/s from harmonic TST.

    See http://dx.doi.org/10.18419/opus-9841, chapter 5. Contrary to
    the Eyring rate this function does only takes a scalar temperature as the
    partition functions are also functions of the temperature and would
    have to be recalculated for different temperatures.

    A possible extension would be to also support multiple rs/ts partition functions,
    one for each temperature.

    Parameters
    ----------
    barrier_height
        Barrier height (energy, enthalpy, gibbs energy, ...) in Hartree.
    rs_part_func
        Partition function of the reactant state.
    ts_part_func
        Partition function of the transition state.
    temperature
        Temperature in Kelvin.
    trans_coeff
        Unitless transmission coefficient, e.g., obtained from Wigner or Eckart
        correction.

    Returns
    -------
    rate
        HTST reaction rate in 1/s.
    """
    rate_eyring = eyring_rate(barrier_height, temperature, trans_coeff)[0]
    rate = ts_part_func / rs_part_func * rate_eyring
    return rate


def wigner_corr(
    temperature: float,
    imag_frequency: float,
) -> float:
    """Tunneling correction according to Wigner.

    See https://doi.org/10.1002/qua.25686 eq. (12) and
    https://doi.org/10.1007/978-3-642-59033-7_9 for the original publication.

    Parameters
    ----------
    temperature
        Temperature in Kelvin.
    imag_frequency
        Imaginary frequency in 1/s.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Wigner.
    """
    kappa = 1 + 1 / 24 * (PLANCK * abs(imag_frequency) / (KB * temperature)) ** 2
    return kappa


def eckart_corr(
    fw_barrier_height: float,
    bw_barrier_height: float,
    temperature: float,
    imag_frequency: float,
) -> float:
    """Tunneling correction according to Eckart.

    See [3], [4] and [5]. The correction should be independent of the order
    fw_barrier_height/bw_barrier_height.

    Parameters
    ----------
    fw_barrier_height
        Barrier height in forward direction in Hartree.
    bw_barrier_height
        Barrier height in backward direction in Hartree.
    temperature
        Temperature in Kelvin.
    imag_frequency
        Frequency in 1/s of the imaginary mode at the TS.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Eckart.
    """
    kBT = KBAU * temperature  # Hartree
    two_pi = 2 * np.pi

    imag_frequency = abs(imag_frequency) * AU2SEC  # Convert from 1/s to 1/au
    hnu = PLANCKAU * imag_frequency
    u = hnu / kBT  # unitless
    v1 = fw_barrier_height / kBT
    v2 = bw_barrier_height / kBT
    alpha1, alpha2 = [
        two_pi * barrier / hnu for barrier in (fw_barrier_height, bw_barrier_height)
    ]
    quot = 1 / (1 / np.sqrt(alpha1) + 1 / np.sqrt(alpha2))
    d = 1 / two_pi * np.sqrt(np.abs(4 * alpha1 * alpha2 - np.pi ** 2))
    cosh_d = np.cosh(two_pi * np.abs(d))

    def eps(E):
        return (E - fw_barrier_height) / kBT

    def ai(E, vi):
        epsilon = eps(E)
        return np.sqrt(2 * (epsilon + vi) / (np.pi * u)) * quot

    def _a1(E):
        return ai(E, vi=v1)

    def _a2(E):
        return ai(E, vi=v2)

    def eckart_probability(E):
        a1 = _a1(E)
        a2 = _a2(E)
        plus = np.cosh(two_pi * (a1 + a2))
        minus = np.cosh(two_pi * (a1 - a2))
        return (plus - minus) / (plus + cosh_d)

    def eckart_func(E):
        epsilon = eps(E)
        P = eckart_probability(E)
        return P * np.exp(-epsilon)

    # Integration limits
    E_low = max(0, fw_barrier_height - bw_barrier_height)
    E_max = 15  # 15 Hartree max
    kappa, _ = integrate.quad(eckart_func, E_low, E_max)
    # Make kappa unitless again by dividing by kBT, as we integrated over the energy.
    kappa /= kBT
    return kappa


def tunl(alph1: float, alph2: float, U: float, strict: bool = False):
    """Eckart correction factor for rate constants according to Brown.

    Python adaption of the TUNL subroutine in 4. Appendix of [5].

    Parameters
    ----------
    alph1
        Unitless barrier height descriptor. 2π V1 / (h nu*); see (2) in [5].
    alph2
        Unitless barrier heigth descriptor. 2π V2 / (h nu*); see (2) in [5].
    u
        Unitless curvature descriptor. h nu* / kT; see (2) in [5].
    strict
        If enabled, arguments are bound checked. Will raise AssertionError if
        they are out of bonds. TUNL was found to yield accurate results when
        the arguments are within bounds.

    Returns
    -------
    G
        Unitless tunneling correction according to Eckart.
    """
    if strict:
        alphs = np.array((alph1, alph2))
        assert (0.5 <= alphs).all() and (alphs <= 20).all() and 0 < U <= 16
        if (alphs >= 8).all():
            assert U <= 12
        if (alphs >= 16).all():
            assert U <= 10

    # Quadrature arguments
    x = np.array((-0.9324695, -0.6612094, -0.2386192, 0.2386192, 0.6612094, 0.9324695))
    w = np.array((0.1713245, 0.3607616, 0.4679139, 0.4679139, 0.3607616, 0.1713245))

    pi = np.pi
    upi2 = U / (2 * pi)

    C = 0.125 * pi * U * (1 / np.sqrt(alph1) + 1.0 / np.sqrt(alph2)) ** 2
    v1 = upi2 * alph1
    v2 = upi2 * alph2
    D = 4 * alph1 * alph2 - pi ** 2
    DF = np.cosh(np.sqrt(D)) if D >= 0 else np.cos(np.sqrt(-D))

    ez = -v1 if (v2 >= v1) else -v2
    eb = min(
        C * (np.log(2.0 * (1.0 + DF) / 0.014) / (2 * pi)) ** 2 - (v1 + v2) / 2, 3.2
    )
    em = (eb - ez) / 2
    ep = (eb + ez) / 2

    # Original loop
    E = em * x + ep
    A1 = pi * np.sqrt((E + v1) / C)
    A2 = pi * np.sqrt((E + v2) / C)
    fm = np.cosh(A1 - A2)
    fp = np.cosh(A1 + A2)
    G = (w * np.exp(-E) * (fp - fm) / (fp + DF)).sum()
    # Outside the loop
    G = em * G + np.exp(-eb)
    return G


def eckart_corr_brown(
    fw_barrier_height: float,
    bw_barrier_height: float,
    temperature: float,
    imag_frequency: float,
) -> float:
    """Tunneling correction according to Eckart.

    Wrapper for the TUNL subroutine as given in the appendix of [5].

    Parameters
    ----------
    fw_barrier_height
        Barrier height in forward direction in Hartree.
    bw_barrier_height
        Barrier height in backward direction in Hartree.
    temperature
        Temperature in Kelvin.
    imag_frequency
        Frequency in 1/s of the imaginary mode at the TS.

    Returns
    -------
    kappa
        Unitless tunneling correction according to Eckart.
    """
    kBT = KBAU * temperature
    hnu = PLANCKAU * imag_frequency * AU2SEC
    u = hnu / kBT

    def alpha(barr):
        return 2 * np.pi * barr / hnu

    alpha1 = alpha(fw_barrier_height)
    alpha2 = alpha(bw_barrier_height)
    kappa = tunl(alpha1, alpha2, u)
    return kappa
