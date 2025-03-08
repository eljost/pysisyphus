import numpy as np

from pysisyphus.constants import KBAU


def boltzmann_weights(energies: np.ndarray, T: float) -> np.ndarray:
    """Boltzmann probability/weights for given energies and temperature.

    Energy in au, temperature in Kelvin."""
    kbT = KBAU * T
    rel_energies = energies - energies.min()
    exp_terms = np.exp(-rel_energies / kbT)
    part_func = np.sum(exp_terms)
    weights = exp_terms / part_func
    return weights
