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
# [4] https://doi.org/10.1103/PhysRevB.62.7809
#     Approximate two-electron spin-orbit coupling term for density-functional-theory
#     DFT calculations using the Douglas-Kroll-Hess transformation
#     Boettger, 2000
# [5] https://doi.org/10.1021/jacsau.2c00659
#     State Interaction Linear Response Time-Dependent
#     Density Functional Theory with Perturbative Spin−Orbit Coupling:
#     Benchmark and Perspectives
#     Liao, Kasper, Jenkins, Yang, Batista, Frisch, Li, 2023
# [6] https://doi.org/10.1063/1.5020079
#     Electron paramagnetic resonance g-tensors from
#     state interaction spin-orbit coupling density matrix
#     renormalization group
#     Sayfutyarova, Chan, 2018
# [7] https://doi.org/10.1063/5.0130868
#     Spin–orbit couplings within spin-conserving and spin-flipping
#     time-dependent density functional theory:
#     Implementation and benchmark calculations
#     Kotaru, Pokhilko, Krylov, 2022


import functools
from typing import Sequence, TYPE_CHECKING

if TYPE_CHECKING:
    from pyscf import gto
import numpy as np
import scipy as sp

from pysisyphus.constants import AU2NU, CAU
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.excited_states import norm_ci_coeffs
from pysisyphus.wavefunction.shells import PySCFShells
from pysisyphus.wavefunction.helpers import BFType, cart_size, sph_size


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

    3d and 4d effective charges are disabled until I figure they are
    only useful for calculation w/ ECPs or w/o them ... this should
    be described in [2].

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
    # elif 21 <= atomic_num <= 30:
        # return (0.385 + 0.025 * (nev - 2)) * atomic_num
    elif atomic_num in (37, 38) or 49 <= atomic_num <= 54:
        return (0.9228 + 0.0017 * nev) * atomic_num
    # 4d Transition metals from [2]
    # elif 39 <= atomic_num <= 48:
        # return (4.680 + 0.060 * (nev - 2)) * atomic_num
    else:
        return float("nan")


def get_original_charge(atomic_num: int) -> float:
    return float(atomic_num)


def Q(l: int) -> int:
    """Maximum number of accumulated electrons for given angular momentum.

    See A006331 in the On-Line Encyclopedia of Integer Sequences."""
    return l * (l + 1) * (2 * l + 1) // 3


@functools.singledispatch
def boettger_factors2d(
    Ls: Sequence[int], Zs: Sequence[int], bf_type: BFType
) -> np.ndarray:
    """Boettger-factors to scale 1-electron SOC-integrals.

    See eq. (18) in [4]."""

    size_func = cart_size if bf_type == BFType.CARTESIAN else sph_size

    factors = list()
    # Calculate one factor for every shell; the factor is repeated depending
    # on the shells' size, so we have one factor per basis function.
    for L, Z in zip(Ls, Zs):
        shell_size = size_func(L)
        factor = Q(L) / Z
        factors.extend([factor] * shell_size)
    factors = np.sqrt(factors)
    factors2d = 1.0 - factors[:, None] * factors[None, :]
    return factors2d


@boettger_factors2d.register
def _(wf: Wavefunction) -> np.ndarray:
    Ls = list()
    Zs = list()
    for shell in wf.shells:
        Ls.append(shell.L)
        Z = shell.atomic_num
        assert Z > 0, (
            "Can't calculate Boettger factor for basis function "
            "that is not centered on an atom!"
        )
        Zs.append(Z)
    return boettger_factors2d(Ls, Zs, wf.bf_type)


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

    # Coupling between singlet ground state and (excited) triplet states.
    # See eq. (9) in [3]. Indices are chosen to be consistent w/ [3].
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
    nsings = len(XpYs) + 1  # Add one as GS is also included
    ntrips = len(XpYt)
    socs = np.concatenate((gs_t, es_t), axis=0).reshape(nsings, ntrips, 3)
    return socs


@singlet_triplet_so_couplings.register
def _(
    wf: Wavefunction, XpYs: np.ndarray, XpYt: np.ndarray, boettger=False
) -> np.ndarray:
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
    charge_func = get_original_charge if boettger else get_effective_charge
    ints_ao = np.zeros((3, nbfs, nbfs))
    # Loop over all atoms, calculate the spin-orbit integrals and accumulate them
    # with the appropriate effective charge into ints_ao
    for i in range(mol.natm):
        mol.set_rinv_origin(mol.atom_coord(i))
        # charge = mol.atom_charge(i)  # Plain charge
        charge = charge_func(mol.atom_charge(i))
        ints_ao += charge * mol.intor(intor_key)

    if boettger:
        bfactors = boettger_factors2d(wf)
        ints_ao *= bfactors

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


def report_so_couplings(socs):
    nsings, ntrips, _ = socs.shape
    # See for example eq. (15) in [7]
    socs2 = np.abs(socs) ** 2
    tot_socs = np.sqrt(socs2.sum(axis=2))

    socs_nu = socs * AU2NU
    tot_socs_nu = tot_socs * AU2NU
    print()

    fmt = "{: >8.2f}"
    table = TablePrinter(
        header=(
            "S",
            "T",
            "Ms-1 re",
            "Ms-1 im",
            "Ms 0 re",
            "Ms 0 im",
            "Ms+1 re",
            "Ms+1 im",
        ),
        col_fmts=("int3", "int3", fmt, fmt, fmt, fmt, fmt, fmt),
    )
    print("                     <S|H_so|T>, complex numbers in cm⁻¹")
    table.print_header()
    for sstate in range(nsings):
        for tstate in range(ntrips):
            sm1, s0, sp1 = socs_nu[sstate, tstate]
            table.print_row(
                (
                    sstate,
                    tstate + 1,
                    sm1.real,
                    sm1.imag,
                    s0.real,
                    s0.imag,
                    sp1.real,
                    sp1.imag,
                )
            )
        table.print_sep()
    print()

    table = TablePrinter(
        header=("S", "T", "Total", "Ms-1", "Ms 0", "Ms+1"),
        col_fmts=(
            "int3",
            "int3",
            "float_short",
            "float_short",
            "float_short",
            "float_short",
        ),
        width=20,
    )
    print(" " * 15 + "sqrt(sum(|<S|H_so|T>|**2))  |<S|H_so|T>| in cm⁻¹")
    table.print_header()
    for sstate in range(nsings):
        for tstate in range(ntrips):
            ts = tot_socs_nu[sstate, tstate]
            sm1, s0, sp1 = np.abs(socs_nu[sstate, tstate])
            table.print_row((sstate, tstate + 1, ts, sm1, s0, sp1))
        table.print_sep()


def get_soc_states(socs, sing_ens, trip_ens):
    # SOC-matrix to be diagonalized. The diagonal contains the singlet
    # energies including the GS. Every triplet state is repeated three times.
    # There are no singlet-singlet and triplet-triplet couplings, so these
    # blocks are just diagonal matrices. The singlet-triplet couplings appear
    # in the order <S|H_SO|T_1,-1>, <S|H_SO|T_1,0>, <S|H_SO|T_1,1>.
    #
    # The resulting matrix is hermitian.

    # ntrips = len(trips_)
    nsings, ntrips, _ = socs.shape
    assert nsings == len(sing_ens)
    assert ntrips == len(trip_ens)
    socs_flat = socs.reshape(nsings, ntrips * 3)
    soc_mat = np.block(
        [
            [np.diag(sing_ens), socs_flat],
            [np.conj(socs_flat).T, np.diag(np.repeat(trip_ens, 3))],
        ]
    )
    w, v = np.linalg.eigh(soc_mat)
    return w, v


def get_soc_tdms(wf, v, socs, XpYs):
    nsings, ntrips, _ = socs.shape
    nstates = nsings + 3 * ntrips
    sing_tdms = wf.get_transition_dipole_moment(XpYs)
    # Singlet-triplet excitations are spin-forbidden
    trip_tdms = np.zeros((ntrips, 3))
    tdms = np.zeros((nstates, 3))
    tdms[1:nsings] = sing_tdms
    tdms[nsings:] = np.repeat(trip_tdms, 3, axis=0)
    soc_tdms = v.T @ tdms
    return sing_tdms, soc_tdms


def get_soc_oscs(exc_ens, soc_tdms):
    assert len(exc_ens) == len(soc_tdms)
    return (
        2.0
        / 3.0
        * exc_ens
        * np.real(np.einsum("ij,ij->i", np.conj(soc_tdms), soc_tdms, optimize="greedy"))
    )


def root_spin_ms_from_row(row, nsings):
    if row < nsings:
        root = row
        spin = 0
        ms = 0
    else:
        row -= nsings
        root = row // 3
        spin = 1
        ms = {0: -1, 1: 0, 2: 1}[row % 3]
    return root, spin, ms


def report_so_states(w, v, nsings, thresh=1e-2):
    w_nu = w * AU2NU
    w_min_nu = w_nu[0]
    w_nu -= w_min_nu
    gs_shift = w_min_nu
    print(highlight_text("Spin-orbit states"))
    print(f"Shift of SO ground state: {gs_shift:.4f} cm⁻¹")
    print()
    print(
        "                   E(cm⁻¹) Weight   Real       Imag       Label Root   Spin  Mₛ"
    )
    for i, eigval in enumerate(w_nu):
        print(f"SOC state {i: >03d}: {eigval: >12.2f}")
        for row, z in enumerate(v[:, i]):
            root, spin, Ms = root_spin_ms_from_row(row, nsings)
            # Report first triplet state as root 1, not root 0
            if spin == 1:
                root += 1
            label = ("S" if spin == 0 else "T") + str(root)

            z = complex(z)
            re = z.real
            im = z.imag
            norm = (re**2 + im**2) ** 0.5
            weight = norm**2
            if norm >= thresh:
                print(
                    f"\t\t\t {weight:.5f} {re: >10.5f} {im: >10.5f} {label: >5s} {root: >5d} {spin: >5d} {Ms: >5d}"
                )


def report_singlet_trans_moms(w, sing_oscs, sing_tdms):
    w0 = w - w.min()
    w0_nu = AU2NU * w0

    print(highlight_text("Singlet-singlet transition moments"))
    fmt = "{: >10.6f}"
    table = TablePrinter(
        header=(
            "From",
            "To",
            "ΔE in cm⁻¹",
            "fosc",
            "μx",
            "μy",
            "μz",
        ),
        col_fmts=(
            "{: >4d}",
            "{: >4d}",
            "{: >12.2f}",
            "{: >7.5f}",
            fmt,
            fmt,
            fmt,
        ),
    )
    table.print_header()
    for i, osc in enumerate(sing_oscs):
        exc_en_nu = w0_nu[i]
        mux, muy, muz = sing_tdms[i]
        table.print_row((0, i + 1, exc_en_nu, osc, mux, muy, muz))


def report_so_trans_moms(w, soc_oscs, soc_tdms):
    w0 = w - w.min()
    w0_nu = AU2NU * w0

    print(highlight_text("Spin-Orbit Transition Moments"))
    fmt = "{: >10.6f}"
    table = TablePrinter(
        header=(
            "From",
            "To",
            "ΔE in cm⁻¹",
            "fosc",
            "μx re",
            "μx im",
            "μy re",
            "μy im",
            "μz re",
            "μz im",
        ),
        col_fmts=(
            "{: >4d}",
            "{: >4d}",
            "{: >12.2f}",
            "{: >7.5f}",
            fmt,
            fmt,
            fmt,
            fmt,
            fmt,
            fmt,
        ),
    )
    table.print_header()
    for i, osc in enumerate(soc_oscs[1:], 1):
        exc_en_nu = w0_nu[i]
        mux, muy, muz = soc_tdms[i]
        table.print_row(
            (
                0,
                i,
                exc_en_nu,
                osc,
                mux.real,
                mux.imag,
                muy.real,
                muy.imag,
                muz.real,
                muz.imag,
            )
        )


def run(wf, Xas, Yas, Xat, Yat, sing_ens, trip_ens, **kwargs):
    # Singlet-singlet excitations
    nsings, *_ = Xas.shape
    Xas, Yas = norm_ci_coeffs(Xas, Yas)
    XpYs = Xas + Yas

    # Singlet-triplet excitations
    Xat, Yat = norm_ci_coeffs(Xat, Yat)
    XpYt = Xat + Yat

    # The way CI coefficients are normalized in pysisyphus mandates multiplication by
    # sqrt(2). This is also discussed in the PySOC paper.
    XpYs = np.sqrt(2) * XpYs
    XpYt = np.sqrt(2) * XpYt

    ################
    # SO-couplings #
    ################

    socs = singlet_triplet_so_couplings(wf, XpYs, XpYt, **kwargs)
    report_so_couplings(socs)
    # The spherical triplet operators can be converted to Cartesian ones via
    #
    # Tx = (socs[..., 0] - socs[..., 2]) / 2.0
    # Ty = (socs[..., 0] + socs[..., 2]) / 2.0j
    # Tz = 1 / np.sqrt(2.0) * socs[..., 1]
    # See eq. (14) to (16) in [6].
    print()

    ##############
    # SOC-States #
    ##############

    w, v = get_soc_states(socs, sing_ens, trip_ens)

    #####################
    # Report SOC states #
    #####################

    report_so_states(w, v, nsings)
    print()

    #############################
    # Transition dipole moments #
    #############################

    Xasn, Yasn = norm_ci_coeffs(Xas, Yas)
    XpYsn = Xasn + Yasn
    sing_tdms, soc_tdms = get_soc_tdms(wf, v, socs, XpYsn)
    # Drop GS energy from sing_ens
    sing_oscs = wf.oscillator_strength(sing_ens[1:], sing_tdms)
    soc_oscs = get_soc_oscs(w, soc_tdms)
    report_singlet_trans_moms(w, sing_oscs, sing_tdms)
    print()
    print(
        "Singlet-Triplet transition dipole moments are zero, as they are "
        "spin-forbidden!"
    )
    print()
    report_so_trans_moms(w, soc_oscs, soc_tdms)
    print()

    return socs, w, v
