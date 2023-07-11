# [1] https://doi.org/10.1080/00268970110089432
#     Distributed multipole analysis Methods and applications
#     Stone, Alderton, 1985
# [2] https://doi.org/10.1016/0009-2614(81)85452-8
#     Distributed multipole analysis, or how to describe a molecular
#     charge distribution
#     Stone, 1981
# [3] https://doi.org/10.1063/1.3430523
#     Convergence of the multipole expansion for 1,2 Coulomb interactions:
#     The modified multipole shifting algorithm
#     Solano, Pendas, Francisco, Blanco, Popelier, 2010
# [4] https://doi.org/10.1063/5.0076630
#     Multi-center decomposition of molecular densities:
#     A mathematical perspective
#     Benda, Cances, Ehrlacher, Stamm, 2022


import dataclasses
import itertools as it
import time
from typing import Tuple

import numpy as np
from scipy.special import binom

from pysisyphus.constants import BOHR2ANG
from pysisyphus.elem_data import nuc_charges_for_atoms
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.ints.multipole3d_sph import multipole3d_sph
from pysisyphus.wavefunction.helpers import lm_iter
from pysisyphus.wavefunction.solid_harmonics import Rlm


_LE_MAX = 2  # Up to quadrupoles
_TWO_LE_MAX = 2 * _LE_MAX
_SQRT2 = 1 / np.sqrt(2.0)
_DIST_THRESH = 1e-6


def get_index_maps(L_max: int = _LE_MAX) -> dict[(int, int), int]:
    """Map between (l, m) keys and 1d indices in the multipole arrays.

    The ordering here must match the ordering in the integral functions.
    """

    index_map = dict()
    index = 0
    for l in range(L_max + 1):
        for m in range(-l, l + 1):
            key = (l, m)
            index_map[key] = index
            index += 1
    return index_map


# Map between 2d (l, m) and 1d-indices.
_INDEX_MAP = get_index_maps()


def get_C(lmax: int = _TWO_LE_MAX) -> dict[(int, int), complex]:
    """Factors for complex->real transformation of multipoles.

    Defined below A9 in [3].
    """
    C = dict()
    for ma in range(-lmax, lmax + 1):
        for mb in range(-lmax, lmax + 1):
            key = (ma, mb)
            if key == (0, 0):
                value = 1.0
            elif abs(ma) != abs(mb):
                value = 0.0
            elif ma > 0 and mb > 0 and ma == mb:
                value = (-1) ** ma * _SQRT2
            elif ma > 0 and (-ma == mb):
                value = _SQRT2
            elif ma < 0 and (-ma == mb):
                value = -1j * (-1) ** mb * _SQRT2
            elif ma < 0 and mb < 0 and ma == mb:
                value = 1j * _SQRT2
            else:
                raise Exception("Safeguard")
            C[key] = value
    return C


C = get_C()
# Complex conjugate values of C
CCONJ = {k: np.conj(v) for k, v in C.items()}


def get_binom_lmkq(L_max: int = _LE_MAX) -> np.ndarray:
    """Bionmial-coefficient prefactor required for translating multipoles.

    See eq. (2.7) in [1] and eq. (11) in [2]. Returns sqrt of the product
    of two binomial coefficients."""
    lmkq = np.zeros((L_max + 1, 2 * L_max + 1, L_max + 1, 2 * L_max + 1))

    for l in range(L_max + 1):
        for m in range(-l, l + 1):
            for k in range(l + 1):
                for q in range(-k, k + 1):
                    binom1 = binom(l + m, k + q)
                    binom2 = binom(l - m, k - q)
                    lmkq[l, m, k, q] = binom1 * binom2
    lmkq = np.sqrt(lmkq)
    return lmkq


def get_W(xyz: np.ndarray, lmkq: np.ndarray, L_max: int = _LE_MAX) -> np.ndarray:
    """Multipole-shift matrix W.

    See eq. (A11) in [3]."""
    num = (np.arange(L_max + 1) * 2 + 1).sum()
    # Complex type is not strictly needed, as the matrix is always real. But adding
    # to a real matrix raises a warning.
    W = np.zeros((num, num), dtype=complex)
    for (l, m), (ldash, mdash) in it.product(lm_iter(L_max), lm_iter(L_max)):
        lm_ind = _INDEX_MAP[(l, m)]
        lmdash_ind = _INDEX_MAP[(ldash, mdash)]
        for m1 in range(-l, l + 1):
            Cmm1 = C[(m, m1)]
            if Cmm1 == 0.0:
                continue
            for m2 in range(-ldash, ldash + 1):
                Cmdashm2 = CCONJ[(mdash, m2)]
                if Cmdashm2 == 0.0:
                    continue
                # prefact = LMKQ[l, m1, ldash, m2]
                prefact = lmkq[l, m1, ldash, m2]
                if prefact == 0.0:
                    continue
                for m3 in range(ldash - l, l - ldash + 1):
                    Cprod = CCONJ[(m3, m1 - m2)]
                    if Cprod == 0.0:
                        continue
                    Cprod *= Cmm1 * Cmdashm2
                    # Index order is consistent with eq. (A11) in [3] (W_{l'm', lm}).
                    W[lmdash_ind, lm_ind] += prefact * Cprod * Rlm(l - ldash, m3, *xyz)
    return np.real(W)


def move_multipoles_faulty(
    L_max: int,
    lmkq: np.ndarray,
    multipoles: np.ndarray,
    P: np.ndarray,
    S: np.ndarray,
) -> np.ndarray:
    """Move multipoles from P to S. THIS FUNCTION IS FAULTY (read below)!

    This function seems to work for charges and dipole moments, but quadrupole
    moments will come out slightly off. I have no idea what is going on ...

    In my opinion, there is a mistake in Eq. (2.7) in [1]. The argument to the
    spherical harmonics should be (P - S) instead of (S - P). Eq. (11) in [2]
    seems to be correct with (S' = S + a) and R_{L-K, M-Q}(a). When S' = P we have
    a = P - S and must pass P - S to R_{L-K, M-Q} instead of S - P.
    """
    moved_multipoles = np.zeros_like(multipoles)
    PS = P - S
    for l in range(L_max + 1):
        for m in range(-l, l + 1):
            moved_key = _INDEX_MAP[l, m]
            for k in range(l + 1):
                for q in range(-k, k + 1):
                    prefact = lmkq[l, m, k, q]
                    if prefact == 0.0:
                        continue
                    key = _INDEX_MAP[k, q]
                    moved_multipoles[moved_key] += (
                        prefact * multipoles[key] * Rlm(l - k, m - q, *PS)
                    )
    return moved_multipoles


def move_multipoles(
    L_max: int,
    lmkq: np.ndarray,
    multipoles: np.ndarray,
    P: np.ndarray,
    S: np.ndarray,
) -> np.ndarray:
    """Move multipoles from P to S using matrix multiplication.

    Follows the approach outlined in [3].
    """
    W = get_W(P - S, lmkq, L_max)
    # Eq. (7) in [3]
    return multipoles @ W


def closest_sites_and_weights(
    org_coords: np.ndarray,
    sites: np.ndarray,
    site_radii: np.ndarray,
    thresh: float = _DIST_THRESH,
):
    """Determine closest expansion sites and associated weights.

    Implements Stone's closest-neigbours-take-it-all-strategy.
    Alternatively, Vigne-Maeder could be implemented (eq. (78) in [4])."""
    # Determine distances to expansion sites
    dists = np.linalg.norm(sites - org_coords[None, :], axis=1)
    # Take possibly different site radii into account
    dists /= site_radii

    # Determine closest site
    closest_index = dists.argmin()
    # Look for other sites that are also/equally close
    delta_dists = dists - dists[closest_index]
    delta_dist_mask = np.abs(delta_dists) <= thresh
    # Indices of close neighbours
    neighbours = np.arange(len(delta_dists))[delta_dist_mask]
    # Distribute the multipoles equally among closest sites
    nneighs = len(neighbours)
    weights = np.full(nneighs, 1 / nneighs)
    return neighbours, weights


def dma(
    wavefunction: Wavefunction,
    sites: np.ndarray,
    site_radii: np.ndarray,
    exp_thresh: float = 37,
    dens_thresh: float = 1e-9,
):
    assert len(sites) == len(site_radii)
    nmultipoles = (2 * np.arange(_LE_MAX + 1) + 1).sum()
    nsites = len(sites)
    # Array that will hold the final distributed multipoles
    distributed_multipoles = np.zeros((nsites, nmultipoles))

    coords3d = wavefunction.coords.reshape(-1, 3)
    # Sum density matrices into a total density
    P = np.sum(wavefunction.P, axis=0)

    shells = wavefunction.shells
    # Prefactors required for redistributing multipoles
    lmkq = get_binom_lmkq()

    reorder_c2s_coeffs = shells.reorder_c2s_coeffs

    # Distribute nuclear charges on to expansion sites
    nuc_charges = nuc_charges_for_atoms(wavefunction.atoms)
    for i, nuc_charge in enumerate(nuc_charges):
        nuc_multipoles = np.zeros(nmultipoles)
        nuc_multipoles[0] = nuc_charge

        nuc_coords = coords3d[i]
        neighbours, weights = closest_sites_and_weights(nuc_coords, sites, site_radii)
        for neigh_index, weight in zip(neighbours, weights):
            S = sites[neigh_index]
            # Don't move if natural- and expansion-center coincide
            if np.linalg.norm(S - nuc_coords) <= _DIST_THRESH:
                moved_multipoles = nuc_multipoles
            else:
                moved_multipoles = move_multipoles(
                    _LE_MAX,
                    lmkq,  # Binomial coefficient prefactors
                    weight * nuc_multipoles,
                    nuc_coords,  # Original/natural center
                    S,  # Target center/DMA site
                )
            # Sum contributions of the different basis function pairs
            distributed_multipoles[neigh_index] += moved_multipoles

    # Loop over unique shell pairs
    for i, shell_a in enumerate(shells):
        La, A, das, axs, a_ind, a_size = shell_a.as_tuple()
        a_slice_cart = slice(a_ind, a_ind + a_size)
        a_slice_sph = shell_a.sph_slice
        a_c2s = reorder_c2s_coeffs[a_slice_sph, a_slice_cart]
        for j, shell_b in enumerate(shells[i:], i):
            Lb, B, dbs, bxs, b_ind, b_size = shell_b.as_tuple()
            AB = A - B
            AB_dist2 = AB.dot(AB)
            b_slice_cart = slice(b_ind, b_ind + b_size)
            b_slice_sph = shell_b.sph_slice
            b_c2s = reorder_c2s_coeffs[b_slice_sph, b_slice_cart]
            # Pick appropriate integral function
            func = multipole3d_sph[(La, Lb)]

            # Slice of the density matrix
            P_slice = P[a_slice_sph, b_slice_sph]
            # Every off-diagonal element of the density matrix contributes two times.
            sym_factor = 2.0 if (i != j) else 1.0

            # Screen out whole shell pairs
            #
            # The default dens_thresh of 1e-9 is very conservative and will
            # probably never trigger for smaller molecules. Anything above 1e-9,
            # e.g. 1e-8 already leads to quadrupole errors for methane (HF/def2-svp)
            # on the order of 1e-8.
            dens_max = max(abs(P_slice.min()), P_slice.max())
            dens_abs = sym_factor * dens_max
            if dens_abs <= dens_thresh:
                continue

            # Loop over primitive pairs.
            for da, ax in zip(das, axs):
                ax_AB_dist = ax * AB_dist2
                for db, bx in zip(dbs, bxs):
                    px = ax + bx
                    # exp-argument of gaussian overlap
                    exp_arg = (bx * ax_AB_dist) / px
                    # When exp-arg is big then 'exp(-exp_arg)' will be very small
                    # and we can skip this primitive pair. The default thresh is 37:
                    # exp(-37) = 8.533e-17, so the screening is very conservative.
                    # Using 'exp_thresh = 10' results in errors up to 1e-6 in the test
                    # cases.
                    if exp_arg >= exp_thresh:
                        continue

                    # Natural center of Gaussian primitive pair
                    P_prim = (ax * A + bx * B) / (ax + bx)

                    # Calculate multipole integrals.
                    # The minus is taken into account here, as well as the symmetry factor.
                    prim_multipoles = -sym_factor * func(ax, da, A, bx, db, B)
                    # Convert to spherical basis functions & reorder, so the basis function
                    # order matches the density matrix.
                    prim_multipoles = np.einsum(
                        "ij,mjk,kl->mil",
                        a_c2s,
                        prim_multipoles,
                        b_c2s.T,
                        optimize="greedy",
                    )

                    # Contract integrals with density matrix elements.
                    contr_multipoles = prim_multipoles * P_slice[None, :, :]
                    contr_multipoles = contr_multipoles.sum(axis=(1, 2))
                    neighbours, weights = closest_sites_and_weights(
                        P_prim, sites, site_radii
                    )
                    for neigh_index, weight in zip(neighbours, weights):
                        S = sites[neigh_index]
                        # Don't move if natural- and expansion-center coincide
                        if np.linalg.norm(S - P_prim) <= _DIST_THRESH:
                            moved_multipoles = contr_multipoles
                        else:
                            moved_multipoles = move_multipoles(
                                _LE_MAX,
                                lmkq,  # Binomial coefficient prefactors
                                weight * contr_multipoles,
                                P_prim,  # Original/natural center
                                S,  # Target center/DMA site
                            )
                        # Add to final distributed multipoles at sites
                        distributed_multipoles[neigh_index] += moved_multipoles
            # End Loop over primitive pairs
    return distributed_multipoles


def get_radii(atoms: Tuple[str]) -> np.ndarray:
    """Default site radii from DMA. 0.325 Å for Hydrogen, 0.65 Å else."""
    radii = [0.325 if atom.lower() == "h" else 0.65 for atom in atoms]
    radii = np.array(radii) / BOHR2ANG
    return radii


def run_dma(wf, sites=None):
    if sites is None:
        dma_sites = wf.coords.reshape(-1, 3).copy()
    else:
        dma_sites = np.array(sites).reshape(-1, 3)
    site_radii = get_radii(wf.atoms)
    if sites is not None:
        assert len(sites) == 3
        site_radii = np.array((1.0,))
    # TODO: handling radius assignment when sites are given ...
    assert len(site_radii) == len(dma_sites)

    start = time.time()
    distributed_multipoles = dma(wf, dma_sites, site_radii)
    dur = time.time() - start
    print(f"Distributed multipole analysis took {dur:.4f} s")
    cfmt = " >12.6f"
    mfmt = " >12.6f"

    for i, dma_site in enumerate(dma_sites):
        x, y, z = dma_site
        radius = site_radii[i]
        # charge = dma_charges[i]
        charge = -999
        print(
            f"(Z={charge}, center: {x=:{cfmt}},  {y=:{cfmt}}, {z=:{cfmt}}, {radius=:{cfmt}}) au"
        )
        site_multipoles = distributed_multipoles[i]
        index = 0
        for l in range(_LE_MAX + 1):
            size = 2 * l + 1
            l_multipoles = site_multipoles[index : index + size]
            l_norm = np.linalg.norm(l_multipoles)
            multi_strs = [
                f"{sm:{mfmt}}" for sm in site_multipoles[index : index + size]
            ]
            print(f"\t|Q_{l}|={l_norm:{mfmt}}: {', '.join(multi_strs)}")
            index += size
        print()
    return distributed_multipoles
