import numpy as np
import numpy.typing as npt

from pysisyphus.constants import KB, KBAU, PLANCK


def eyring_rate(
    barrier_height: float,
    temperature: npt.ArrayLike,
    trans_coeff: float = 1.0,
) -> npt.ArrayLike:
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
        Transmission coefficient, unitless.

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

    See http://dx.doi.org/10.18419/opus-9841, chapter 5.

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
        Transmission coefficient, unitless.

    Returns
    -------
    rate
        HTST reaction rate in 1/s.
    """
    rate_eyring = eyring_rate(barrier_height, temperature, trans_coeff)
    rate = ts_part_func / rs_part_func * rate_eyring
    return rate
