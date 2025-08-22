# [1] https://doi.org/10.1063/5.0021923
#     Extending nudged elastic band method to reaction pathways
#     involving multiple spin states
#     Zhao, Watanabe, Nakatani, Nakayama, Xu, Hasegawa, 2020

import numpy as np


def middle_is_meisc(energies: np.ndarray, allow_fast_flip: bool = False) -> bool:
    """Determine if the middle pair of three energy pairs is an MEISC candidate.

    See Discussion in Section IV.C of [1].

    Paramters
    ---------
    energies
        2d array of shape (3, 2), containing the energies of 2 states for 3 images.
    allow_fast_flip
        Boolean flag, controlling whether the spin can flip two times in triple.

    Returns
    -------
    is_meisc
        Boolean that indicates whether the middle image is an MEISC candidate.
    """
    assert energies.shape == (3, 2)
    # Determine indices (0 or 1) of the states with lowest energy.
    states = list(np.argmin(energies, axis=1))
    gaps = np.abs(energies[:, 0] - energies[:, 1])

    # See fig. [6] in [1]
    match states:
        # Case a
        case [0, 1, 1] | [1, 1, 0]:
            is_meisc = gaps[1] < gaps[0]
        # Case b
        case [0, 0, 1] | [1, 0, 0]:
            is_meisc = gaps[1] < gaps[2]
        # By default we assume that a spin flip can't happen so fast
        # (allow_fast_flip = False). This may not be justified for widely
        # spaced images.
        #
        # Case c in the paper.
        case [0, 1, 0] | [1, 0, 1]:
            is_meisc = allow_fast_flip
        # Default; fallback
        case _:
            is_meisc = False
    return is_meisc


def determine_meisc_images(all_energies: np.ndarray, **kwargs) -> list[int]:
    """Determine MEISC image indices for an array of energy pairs.

    Parameters
    ----------
    all_energies
        2d array of shape (nimages, nstates). nstates must be 2.
    **kwargs
        Additional kwargs that will be passed to 'middle_is_meisc()'

    Returns
    -------
    meisc_inds
        List of integer images of MEISC images.
    """
    _, nstates = all_energies.shape
    assert nstates == 2
    meisc_inds = []
    nimages = len(all_energies)
    i = 0
    while True:
        if i >= nimages:
            break
        triples = all_energies[i : i + 3]
        # Break when not enough images are left for a full triple
        if len(triples) < 3:
            break
        # If we obtain an MEISC index we skip triple generation starting from the
        # immediate neighbour.
        if middle_is_meisc(triples, **kwargs):
            meisc_inds.append(i + 1)
            inc = 2
        # If no MEISC index is obtained we also generate a triple starting from the
        # immediate neighbour.
        else:
            inc = 1
        i += inc
    return meisc_inds
