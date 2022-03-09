import numpy as np
import numpy.typing as npt

from pysisyphus.constants import KB, KBAU, PLANCK


def eyring_rate(
    barrier_height: npt.ArrayLike,
    temperature: float,
    trans_coeff: float = 1.0,
) -> npt.ArrayLike:
    """Rate constant in 1/s from the Eyring equation.

    See https://pubs.acs.org/doi/10.1021/acs.organomet.8b00456
    "Reaction Rates and Barriers" on p. 3234 and eq. (8).

    Parameters
    ----------
    barrier_height : float
        Barrier height (energy, enthalpy, gibbs energy, ...) in Hartree.
    temperature : float
        Temperature in Kelvin.
    trans_coeff : float, optional, default = 1.0
        Transmission coefficient.

    Returns
    -------
    rate : float
        Reaction rate in 1/s.
    """
    barrier_height = np.array(barrier_height, dtype=float)
    prefactor = trans_coeff * KB * temperature / PLANCK
    rate = prefactor * np.exp(-barrier_height / (KBAU * temperature))
    return rate
