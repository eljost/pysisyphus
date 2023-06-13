# [1] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.7b00325
#     Franck−Condon Models for Simulating the Band Shape of
#     Electronic Absorption Spectra
#     Li, Truhlar, 2017


import numpy as np


def lq2_abs_cross_sec(
    E: np.ndarray,
    f_osc: float,
    dE_vert: float,
    angfreqs: np.ndarray,
    displacements: np.ndarray,
):
    """LQ2 absorption cross section.

    Eq. (11) in [1].

    Parameters
    ----------
    E
        Array containing the energies of the incident photon in atomic units.
    f_osc
        Oscillator strength.
    dE_vert
        Vertic excitation energy of the excited state in atomic units.
    angfreqs
        Array of normal mode angular frequncies in atomic units.
    displacements
        Array of unitless displacements in along the normal modes.

    Returns
    -------
    sigma
        Array containing absorption cross sections. Its maximum will be
        centered at the excitation energy. In [1] the peak is shifted by
        its HWHM. See just below eq. (11) in [1] for the formula. While
        the sigmas are calculated for the values in the array 'E', they
        should be plotted at 'E' + 'hwhm'.
    """
    assert len(angfreqs) == len(displacements)

    denom = (displacements**2 * angfreqs**2).sum()
    num = (E - dE_vert) ** 2
    sigma = f_osc * np.exp(-num / denom)
    return sigma
