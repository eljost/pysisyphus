# [1] https://doi.org/10.1002/jcc.21113
#     One-electron spin-orbit contribution by effective nuclear charges
#     Chiodo, Russo, 2008
# [2] https://doi.org/10.1021/jp983453n
#     Effective Nuclear Charges for the First- through Third-Row
#     Transition Metal Elements in Spin−Orbit Calculations
#     Koseki, Schmidt, Gordon, 1998
# [3] https://doi.org/10.1021/acs.jctc.6b00915
#     Evaluation of Spin-Orbit Couplings with Linear-Response
#     Time-Dependent Density Functional Methods
#     Gao, Bai, Fazzi, Niehaus, Barbatti, Thiel, 2017

import functools
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyscf import gto
import numpy as np
import scipy as sp

from pysisyphus.constants import CAU
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.shells import PySCFShells
from pysisyphus.wavefunction.helpers import BFType


# Prefactor in atomic units in eq. (1) in [3]; e and m_e are 1.0
_FACTOR = 1 / (2 * CAU**2)


def get_cart_norms(mol: "gto.Mole") -> np.ndarray:
    """Get normalization matrix to ensure self-overlaps of 1.0.

    In PySCF not all Cartesian basis functions have a self-overlap of 1.0.
    This can be fixed by this matrix.

    Parameters
    ----------
    mol
        pyscf.gto.Mole object.

    Returns
    -------
    NN
        2d matrix of shape (nao, nao) containing normalization factors
        for 2-center integrals over Cartesian basis functions.
    """
    S = mol.intor("int1e_ovlp_cart")
    N = 1 / np.diag(S) ** 0.5
    NN = N[:, None] * N[None, :]
    return NN


def get_pyscf_P_sph(shells):
    sph_Ps = PySCFShells.sph_Ps
    P_sph = sp.linalg.block_diag(*[sph_Ps[shell.L] for shell in shells])
    return P_sph


def get_effective_charge(atomic_num: int) -> float:
    """Effective charge for SOC calculations as described in [1].

    Parameter
    ---------
    atomic_num
        Atomic number, positive integer.

    Returns
    -------
    Zeff
        Effective charge.
    """

    # Number of valence electrons for given atomic number
    # fmt: off
    ne_valence = [
        1, 2,  # 1st period
        1, 2, 3, 4, 5, 6, 7, 8,  # 2nd period
        1, 2, 3, 4, 5, 6, 7, 8,  # 3rd period
        1, 2,  # 4th period, main group
        3, 4, 5, 6, 7, 8, 9, 10, 11, 12,  # 3d-TM
        3, 4, 5, 6, 7, 8,  # 4th period, main group
        1, 2,  # 5th period, main group
        3, 4, 5, 6, 7, 8, 9, 10, 11, 12,  # 4d-TM
        3, 4, 5, 6, 7, 8,  # 5th period, main group
    ]
    # As indices start at 0 but atomic numbers at 1 we substract 1.
    nev = ne_valence[atomic_num - 1]

    if atomic_num == 1:
        return 1.0
    elif atomic_num == 2:
        return 2.0
    elif 3 <= atomic_num <= 10:
        return (0.2517 + 0.0626 * nev) * atomic_num
    elif 10 <= atomic_num <= 18:
        return (0.7213 + 0.0144 * nev) * atomic_num
    elif atomic_num in (19, 20) or 31 <= atomic_num <= 36:
        return (0.8791 + 0.0039 * nev) * atomic_num
    # 3d Transition metals from [2]
    elif 21 <= atomic_num <= 30:
        return (0.385 + 0.025 * (nev - 2)) * atomic_num
    elif atomic_num in (37, 38) or 49 <= atomic_num <= 54:
        return (0.9228 + 0.0017 * nev) * atomic_num
    # 4d Transition metals from [2]
    elif 39 <= atomic_num <= 48:
        return (4.680 + 0.060 * (nev - 2)) * atomic_num
    else:
        return float("nan")


@functools.singledispatch
def singlet_triplet_so_couplings(
    C: np.ndarray, ints_ao: np.ndarray, XpYs: np.ndarray, XpYt: np.ndarray
) -> np.ndarray:
    """Singlet-triplet spin-orbit couplings from LR-TDDFT.

    As desribed in [3].

    Parameters
    ----------
    C
        MO-coefficients, 2d array with shape (nao, nmo).
    ints_ao
        Spin-orbit integrals in AO basis with shape (3, nao, nao).
    XpYs
        X and Y vectors for singlet-singlet exictations from TD-DFT.
        3d array with shape (nsings, nocc, nvirt).
    XpYs
        X and Y vectors for singlet-triplet exictations from TD-DFT.
        3d array with shape (ntrips, nocc, nvirt).

    Returns
    -------
    socs
        2d array with shape ((nsings + 1) * ntrips, 3) containing the
        spin-orbit couplings in atomic units.
    """
    nao, _ = C.shape
    assert ints_ao.shape == (3, nao, nao)
    assert XpYs.shape[1:] == XpYt.shape[1:]

    # Transform AO integrals to MO basis
    ints_mo = C.T @ ints_ao @ C
    ints_mo = _FACTOR * ints_mo

    # Determine number of active occupied and virtual orbitals from the shape of the
    # CI-coefficient matrix.
    _, occ_act, virt_act = XpYs.shape

    # Singlet ground state to triplet states, eq. (9) w/ consistent indices
    ints_act = ints_mo[:, :occ_act, occ_act : occ_act + virt_act]
    Ch_gs = np.einsum("Ijb,kjb->kI", XpYt, ints_act, optimize="greedy")
    # Coupling Singlet - T_1,1
    gs_tp1 = -1 / 2 * (Ch_gs[0] + (1j * Ch_gs[1]))
    # Coupling Singlet - T_1,-1
    gs_tm1 = -gs_tp1
    # Coupling GS Singlet - T_1,0
    gs_t0 = 1 / np.sqrt(2) * Ch_gs[2]
    gs_t = np.stack((gs_tm1, gs_t0, gs_tp1), axis=1)

    # Excited singlet states to triplet states Eq. (8)
    ints_act_occ = ints_mo[:, :occ_act, :occ_act]
    # In the PySOC code the case i == j is skipped
    dr, dc = np.diag_indices(occ_act)
    ints_act_occ[:, dr, dc] = 0.0
    Ch_es_occ = np.einsum(
        "Iia,Jja,kji->kIJ", XpYs, XpYt, ints_act_occ, optimize="greedy"
    )
    ints_act_virt = ints_mo[:, occ_act:, occ_act:]
    # In the PySOC code the case a == b is skipped
    dr, dc = np.diag_indices(virt_act)
    ints_act_virt[:, dr, dc] = 0.0
    Ch_es_virt = np.einsum(
        "Iia,Jib,kab->kIJ", XpYs, XpYt, ints_act_virt, optimize="greedy"
    )
    # Coupling Singlet - T_1,1
    es_tp1 = (
        1.0
        / (2 * np.sqrt(2.0))
        * ((Ch_es_occ[0] + 1j * Ch_es_occ[1]) - (Ch_es_virt[0] + 1j * Ch_es_virt[1]))
    )
    # Coupling Singlet - T_1,-1
    es_tm1 = -es_tp1
    # Coupling Singlet - T_1,0
    es_t0 = 1.0 / 2.0 * (-Ch_es_occ[2] + Ch_es_virt[2])
    es_t = np.stack((es_tm1.flatten(), es_t0.flatten(), es_tp1.flatten()), axis=1)

    # Fuse GS-ES and ES-ES spin-orbit couplings
    socs = np.concatenate((gs_t, es_t), axis=0)
    return socs


@singlet_triplet_so_couplings.register
def _(wf: Wavefunction, XpYs: np.ndarray, XpYt: np.ndarray) -> np.ndarray:
    """Wrapper that prepares all required quantites from pysisyphus WF and PySCF."""

    assert wf.restricted, "Unrestricted calculations are currently not supported!"

    shells = wf.shells
    cartesian = wf.bf_type == BFType.CARTESIAN
    if cartesian:
        nbfs = shells.cart_size
        intor_key = "int1e_prinvxp_cart"
        # PySCF has a sensible ordering of Cartesian basis functions.
        # xx, xy, xz, yy, yz, zz ... lexicographic.
        perm_pyscf = np.eye(nbfs)
    else:
        nbfs = shells.sph_size
        intor_key = "int1e_prinvxp_sph"
        perm_pyscf = get_pyscf_P_sph(shells)

    # Permutation matrix to reorder the MO coefficients
    perm = wf.get_permut_matrix()

    # Create PySCF mol from pysisyphus shells object
    mol = shells.to_pyscf_mol()

    # Calculate spin-orbit integrals w/ PySCF
    ints_ao = np.zeros((3, nbfs, nbfs))
    # Loop over all atoms, calculate the spin-orbit integrals and accumulate them
    # with the appropriate effective charge into ints_ao
    for i in range(mol.natm):
        mol.set_rinv_origin(mol.atom_coord(i))
        # charge = mol.atom_charge(i)  # Plain charge
        charge = get_effective_charge(mol.atom_charge(i))
        ints_ao += charge * mol.intor(intor_key)

    # Normalize Cartesian bfs with L >= 2, as they don't have unit self-overlaps in PySCF.
    if cartesian:
        N = get_cart_norms(mol)
        ints_ao = N * ints_ao

    # Bring AO integrals from PySCF into pysisphus-order.
    # For Cartesian basis functions the permutation matrix is the unit matrix,
    # but for spherical basis functions this will have an effect on the basis
    # function order.
    ints_ao = perm_pyscf.T @ ints_ao @ perm_pyscf

    # Considering alpha MO coefficients is enough as we have a restricted wavefunction.
    Ca, _ = wf.C
    # Reorder MO-coefficients from external order to pysisyphus order.
    Ca = perm.T @ Ca

    return singlet_triplet_so_couplings(Ca, ints_ao, XpYs, XpYt)
