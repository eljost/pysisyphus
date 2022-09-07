from typing import Iterable, List, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray


def mo_inds_from_tden(
    T: NDArray[float], ci_thresh: float
) -> Tuple[NDArray[int], NDArray[int], NDArray[float]]:
    """Indices (from, to) for CI-coefficients above a threshold."""
    assert T.ndim == 2
    occ, _ = T.shape
    from_, to_ = np.where(np.abs(T) > ci_thresh)
    coeffs = T[from_, to_]
    to_ += occ
    return from_, to_, coeffs


def to_pairs(arr_a: Iterable, arr_b: Iterable) -> List[frozenset]:
    """Zip two iterables to a list of frozensets."""
    return [frozenset(sorted(p)) for p in zip(arr_a, arr_b)]


def to_same_basis(
    pairs1: Sequence[frozenset],
    vals1: Iterable,
    pairs2: Sequence[frozenset],
    vals2: Iterable,
) -> Tuple[List[int], List[int], NDArray, NDArray]:
    unique_pairs = tuple(set(pairs1 + pairs2))
    zero_vec = np.zeros(len(unique_pairs))
    pair_map = {p: i for i, p in enumerate(unique_pairs)}

    def sort(pairs, vals):
        vals_sorted = zero_vec.copy()
        for pair, v in zip(pairs, vals):
            vals_sorted[pair_map[pair]] = v
        return vals_sorted

    vals1_sorted = sort(pairs1, vals1)
    vals2_sorted = sort(pairs2, vals2)
    try:
        from_, to_ = zip(*unique_pairs)
    # Raised when no unique pairs are present
    except ValueError:
        from_, to_ = list(), list()
    return from_, to_, vals1_sorted, vals2_sorted


def ovlp(
    from_A: NDArray[int],
    to_A: NDArray[int],
    vecA: NDArray[float],
    from_B: NDArray[int],
    to_B: NDArray[int],
    vecB: NDArray[float],
    S_MO: NDArray[float],
) -> NDArray[float]:
    """Overlap between two transition orbital pairs.

    According to Eq. (9)"""
    numA = vecA.size
    numB = vecB.size
    S = np.zeros((numA, numB))

    for i, coeff_a in enumerate(vecA):
        for j, coeff_b in enumerate(vecB):
            s_occ = S_MO[from_A[i], from_B[j]]
            s_vir = S_MO[to_A[i], to_B[j]]
            # Eq. (9)
            S[i, j] = coeff_a / abs(coeff_a) * coeff_b / abs(coeff_b) * s_occ * s_vir
    return S


def r_diff(
    vecP_from: NDArray[float], vec_to: NDArray[float], t_to: NDArray[float]
) -> float:
    abs_vec_to = np.abs(vec_to)

    def diff(vec_to):
        """Eq. (14)."""
        return np.abs(np.abs(vecP_from + vec_to) * t_to).sum()

    r_plus = diff(abs_vec_to)
    r_minus = diff(-1 * abs_vec_to)
    r = min(r_plus, r_minus)
    return r


def parse_tden(
    T: NDArray[float], thresh: float
) -> Tuple[List[frozenset], NDArray[int], NDArray[int], NDArray[float]]:
    from_T, to_T, coeffs = mo_inds_from_tden(T, thresh)
    pairs = to_pairs(from_T, to_T)
    return pairs, from_T, to_T, coeffs


def rAB(
    TXA: NDArray[float],
    TYA: NDArray[float],
    TXB: NDArray[float],
    TYB: NDArray[float],
    S_MO: NDArray[float],
    exc_thresh: float = 0.05,
    deexc_thresh: float = 0.02,
) -> float:
    """Transition-orbital-pair overlaps.

    See 'Transition orbital projection approach for excited state tracking' in
    J. Chem. Phys. 156, 214104 (2022); https://doi.org/10.1063/5.0081207

    Parameters
    ----------
    TXA
        Excitation coefficient (X) vector for 1 state at geometry A.
    TYA
        Deexcitation coefficient (Y) vector for 1 state at geometry A.
    TXB
        Excitation coefficient (X) vector for 1 state at geometry B.
    TYB
        Deexcitation coefficient (Y) vector for 1 state at geometry B.
    S_MO
        Matrix containing overlaps between the MOs at geometry A and B.
    exc_thresh
        Coefficient threshold up to which CI-coefficients from TXA and TXA
        are still considered.
    deexc_thresh
        Coefficient threshold up to which CI-coefficients from TYA and TYA
        are still considered.
    """
    # A
    a_pairs_A, a_from_A, a_to_A, alphaA = parse_tden(TXA, exc_thresh)
    b_pairs_A, b_from_A, b_to_A, betaA = parse_tden(TYA, deexc_thresh)
    a_from_A, a_to_A, alphaA, betaA = to_same_basis(a_pairs_A, alphaA, b_pairs_A, betaA)
    b_from_A = a_from_A
    b_to_A = a_to_A
    # B
    a_pairs_B, a_from_B, a_to_B, alphaB = parse_tden(TXB, exc_thresh)
    b_pairs_B, b_from_B, b_to_B, betaB = parse_tden(TYB, deexc_thresh)
    a_from_B, a_to_B, alphaB, betaB = to_same_basis(a_pairs_B, alphaB, b_pairs_B, betaB)
    b_from_B = a_from_B
    b_to_B = a_to_B

    # Eqs. (7) and (8)
    XA = alphaA + betaA
    YA = alphaA - betaA
    XB = alphaB + betaB
    YB = alphaB - betaB

    # Eqs. (10) and (11)
    SX = ovlp(a_from_A, a_to_A, XA, a_from_B, a_to_B, XB, S_MO)
    XPA = (np.abs(XA)[None, :] @ SX).flatten()
    XPB = (SX @ np.abs(XB)[:, None]).flatten()

    SY = ovlp(b_from_A, b_to_A, YA, b_from_B, b_to_B, YB, S_MO)
    YPA = (np.abs(YA)[None, :] @ SY).flatten()
    YPB = (SY @ np.abs(YB)[:, None]).flatten()

    # Eq. (13); Where does the 2 come from?!
    tA = 2 * (alphaA ** 2 - betaA ** 2)
    tB = 2 * (alphaB ** 2 - betaB ** 2)

    rXAB = r_diff(XPA, XB, tB)  # r_X^A->B
    rYAB = r_diff(YPA, YB, tB)  # r_Y^A->B
    rXBA = r_diff(XPB, XA, tA)  # r_X^B->A
    rYBA = r_diff(YPB, YA, tA)  # r_Y^B->A
    # Eq. (15)
    r = 0.5 * (rXAB + rYAB + rXBA + rYBA)
    return r
