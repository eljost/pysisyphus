# [1] https://doi.org/10.1016/j.ccr.2018.01.019
#     Quantitative wave function analysis for excited states
#     of transition metal complexes
#     Mai et al., 2018

import itertools as it
from functools import partial
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from pysisyphus.wavefunction.wavefunction import Wavefunction


def norm_ci_coeffs(
    Xa: NDArray[float],
    Ya: NDArray[float],
    Xb: Optional[NDArray[float]] = None,
    Yb: Optional[NDArray[float]] = None,
    restricted_norm: float = 0.5,
    unrestricted_norm: float = 1.0,
) -> Tuple[NDArray[float], NDArray[float]]:
    """Normalize transition density matrices.

    target_norm = (N*X)**2 - (N*Y)**2
                = N**2 * (X**2 - Y**2)
    N**2 = target_norm / (X**2 - Y**2)
    N    = sqrt(target_norm / (X**2 - Y**2))
    """
    nstates_a, occ_a, virt_a = Xa.shape
    unrestricted = (Xb is not None) and (Yb is not None)
    if unrestricted:
        nstates_b, occ_b, virt_b = Xb.shape
        assert nstates_a == nstates_b
        assert (occ_a + virt_a) == (occ_b + virt_b)
        nstates = nstates_a
        X = np.concatenate((Xa.reshape(nstates, -1), Xb.reshape(nstates, -1)), axis=1)
        Y = np.concatenate((Ya.reshape(nstates, -1), Yb.reshape(nstates, -1)), axis=1)
        target_norm = unrestricted_norm
    else:
        nstates = nstates_a
        X = Xa.reshape(nstates, -1)
        Y = Ya.reshape(nstates, -1)
        target_norm = restricted_norm
    X_norms, Y_norms = [np.linalg.norm(mat, axis=1) for mat in (X, Y)]
    ci_norms = X_norms**2 - Y_norms**2
    N = np.sqrt(target_norm / ci_norms)
    X *= N[:, None]
    Y *= N[:, None]
    entries_a = occ_a * virt_a
    Xa = X[:, :entries_a].reshape(nstates, occ_a, virt_a)
    Ya = Y[:, :entries_a].reshape(nstates, occ_a, virt_a)
    returns = [Xa, Ya]
    if unrestricted:
        Xb = X[:, entries_a:].reshape(nstates, occ_b, virt_b)
        Yb = Y[:, entries_a:].reshape(nstates, occ_b, virt_b)
        returns.extend([Xb, Yb])
    return returns


def get_state_to_state_transition_density(
    T_a: NDArray[float], T_b: NDArray[float]
) -> NDArray[float]:
    """State-to-state transition density.

    Parameters
    ----------
    T_a
        Transition density matrix for state A of shape (occ, virt).
    T_b
        Transition density matrix for state B of shape (occ, virt).

    Returns
    -------
    T_ab
        State-to-state transition density of shape (occ+virt, occ+virt).
    """
    assert T_a.shape == T_b.shape

    electron = T_a.T @ T_b
    hole = T_b @ T_a.T
    occ, virt = T_a.shape
    ov = occ + virt
    T_ab = np.zeros((ov, ov))
    T_ab[:occ, :occ] = -hole
    T_ab[occ:, occ:] = electron
    return T_ab


def ct_numbers(
    C: NDArray[float], occ: int, Ts: NDArray[float], nfrags: int, frag_ao_map: Dict
) -> NDArray[float]:
    frag_inds = list(range(nfrags))
    T_full = np.zeros_like(C)
    u, _, vh = np.linalg.svd(C)
    lowdin = u @ vh  # (PQ^T) in eq. (17) of [1]

    omegas = list()
    for T_mo in Ts:
        T = T_full.copy()
        T[:occ, occ:] = T_mo
        Torth = lowdin @ T @ lowdin.T  # eq. (17) in [1]
        Torth2 = Torth**2
        state_omegas = list()
        # Pairs of fragments
        for i, j in it.product(frag_inds, frag_inds):
            aos_i = frag_ao_map[i]
            aos_j = frag_ao_map[j]
            om_ij = Torth2[aos_i][:, aos_j].sum()  # eq. (16) in [1]
            state_omegas.append(om_ij)
        omegas.append(state_omegas)
    return np.array(omegas)


def ct_numbers_for_states(
    wf: Wavefunction,
    frags: List[List[int]],
    Ta: NDArray[float],
    Tb: Optional[NDArray[float]] = None,
):
    """
    Charge-Transfer number for a given wavefunction and fragments.

    Parameters
    ----------
    wf
        Wavefunction object.
    frags
        List of lists of integers, containing the atom indices of the respective
        fragments.
    Ta
        Alpha-spin transition density matrices of shape (nstates, occ_a, virt_a).
    Tb
        Beta-spin transition density matrices of shape (nstates, occ_a, virt_a).

    Returns
    -------
    omegas
        Array of shape (nstates, len(frags)**2), containing the charge-transfer
        numbers of every state.
    """
    shells = wf.shells
    frag_ao_map = shells.fragment_ao_map(frags)
    nfrags = len(frags)

    occ_a, occ_b = wf.occ
    Ca, Cb = wf.C

    ctn = partial(ct_numbers, nfrags=nfrags, frag_ao_map=frag_ao_map)
    om_a = ctn(Ca, occ_a, Ta)
    if Tb is None:
        om_b = om_a
    else:
        om_b = ctn(Cb, occ_b, Tb)
    omegas = om_a + om_b
    return omegas


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

    # Eq. (13) in the original paper includes a factor 2, probably because the algorithm
    # is only formulated for closed-shell cases. Here, we exclude this factor and expect
    # that closed-shell/open-shell transition density matrices are normalized differently,
    # e.g. as done in 'norm_ci_coeffs'
    tA = alphaA**2 - betaA**2
    tB = alphaB**2 - betaB**2

    rXAB = r_diff(XPA, XB, tB)  # r_X^A->B
    rYAB = r_diff(YPA, YB, tB)  # r_Y^A->B
    rXBA = r_diff(XPB, XA, tA)  # r_X^B->A
    rYBA = r_diff(YPB, YA, tA)  # r_Y^B->A
    # Eq. (15)
    r = 0.5 * (rXAB + rYAB + rXBA + rYBA)
    return r


def top_differences(XA, YA, XB, YB, S_MO):
    states_A = XA.shape[0]
    states_B = XB.shape[0]

    rs = list()
    for xai, yai in zip(XA, YA):
        for xbj, ybj in zip(XB, YB):
            r = rAB(xai, yai, xbj, ybj, S_MO)
            rs.append(r)
    rs = np.array(rs).reshape(states_A, states_B)
    return rs


def tden_overlaps(
    C_1: NDArray[float],
    ci_coeffs1: NDArray[float],
    C_2: NDArray[:float],
    ci_coeffs2: NDArray[float],
    S_AO: NDArray[float],
):
    """
    Transition density overlaps.

    Parameters
    ----------
    mo_coeffs1 : ndarray, shape (MOs, AOs)
        MO coefficient matrix. MOs are expected to be in columns.
    mo_coeffs2 : ndarray
        See mo_coeffs1.
    ci_coeffs1 : ndarray, shape(occ. MOs, MOs)
        CI-coefficient matrix.
    ci_coeffs2 : ndarray, shape(occ. MOs, MOs)
        See ci_coeffs1.
    S_AO : ndarray, shape(AOs1, AOs2)
        Double molcule AO overlaps.
    """
    # Total number of molecular orbitals for ci_coeffs1 and ci_coeffs2 (occ + virt)
    nmo1, nmo2 = [sum(ci_coeffs.shape[1:]) for ci_coeffs in (ci_coeffs1, ci_coeffs2)]
    assert S_AO.shape == (nmo1, nmo2)
    _, occ1, _ = ci_coeffs1.shape
    _, occ2, _ = ci_coeffs2.shape
    # MO overlaps, and the respective sub-matrices (occ x occ), (virt x virt)
    S_MO = C_1.T @ S_AO @ C_2
    S_MO_occ = S_MO[:occ1, :occ2]
    S_MO_vir = S_MO[occ1:, occ2:]

    # Thanks Philipp and Klaus!
    overlaps = np.zeros((ci_coeffs1.shape[0], ci_coeffs2.shape[0]))
    for i, state1 in enumerate(ci_coeffs1):
        precontr = S_MO_vir.T @ state1.T @ S_MO_occ
        for j, state2 in enumerate(ci_coeffs2):
            overlaps[i, j] = np.trace(precontr @ state2)

    """
    overlaps = np.einsum(
        "mil,ij,njk,kl->mn",
        ci_coeffs1,
        S_MO_occ,
        ci_coeffs2,
        S_MO_vir.T,
        optimize=["einsum_path", (0, 3), (1, 2), (0, 1)],
    )
    """
    return overlaps


###############################
# Natural Transition Orbitals #
###############################


"""
def nto_overlaps(ntos_1, ntos_2, ao_ovlp):
    states1 = len(ntos_1)
    states2 = len(ntos_2)
    ovlps = np.zeros((states1, states2))
    for i in range(states1):
        n_i = ntos_1[i]
        l_i = n_i.lambdas[:, None]
        ntos_i = l_i * n_i.ntos
        for j in range(states2):
            n_j = ntos_2[j]
            l_j = n_j.lambdas[:, None]
            ntos_j = l_j * n_j.ntos
            ovlp = np.sum(np.abs(ntos_i.dot(ao_ovlp).dot(ntos_j.T)))
            ovlps[i, j] = ovlp
    return ovlps

def nto_org_overlaps(ntos_1, ntos_2, ao_ovlp, nto_thresh=0.3):
    states_1 = len(ntos_1)
    states_2 = len(ntos_2)
    ovlps = np.zeros((states_1, states_2))

    for i in range(states_1):
        n_i = ntos_1[i]
        l_i = n_i.lambdas[:, None]
        ntos_i = n_i.ntos[(l_i >= nto_thresh).flatten()]
        l_i_big = l_i[l_i >= nto_thresh]
        for j in range(states_2):
            n_j = ntos_2[j]
            l_j = n_j.lambdas[:, None]
            ntos_j = n_j.ntos[(l_j >= nto_thresh).flatten()]
            ovlp = np.sum(
                l_i_big[:, None] * np.abs(ntos_i.dot(ao_ovlp).dot(ntos_j.T))
            )
            ovlps[i, j] = ovlp
    return ovlps
"""
