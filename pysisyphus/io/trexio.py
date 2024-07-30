from pathlib import Path
from typing import Optional

import numpy as np
import trexio

from pysisyphus.wavefunction.helpers import BFType
from pysisyphus.wavefunction.normalization import factorial2
from pysisyphus.wavefunction.shells import Shell
from pysisyphus.wavefunction.shells_trexio import TrexIOShells
from pysisyphus.wavefunction import Wavefunction


def norm_pgto(exponent: float, L: int):
    """Primitive normalization factor depending only on L = l + m + n.

    TREXIO seems to use only one normalization factor for every primitive,
    even for L > 0, where one primitive contributes to multiple basis functions.
    """
    L_2 = L / 2
    norm = (
        2**L_2
        * (1 / (2 * exponent)) ** (-L_2 - 0.75)
        / (np.pi**0.75)
        / np.sqrt(factorial2(2 * L - 1))
    )
    return norm


def shells_from_trexio(
    fn: Optional[str | Path] = None, handle: Optional[trexio.File] = None, **kwargs
) -> TrexIOShells:
    assert (fn is not None) or (
        handle is not None
    ), "Either 'fn' or 'handle' must be provided!"

    if fn is not None:
        tf = trexio.File(str(fn), mode="r", back_end=trexio.TREXIO_HDF5)
    else:
        tf = handle
    shells = None

    basis_type = trexio.read_basis_type(tf)
    assert basis_type == "Gaussian"

    shell_num = trexio.read_basis_shell_num(tf)
    nucleus_index = trexio.read_basis_nucleus_index(tf)
    shell_ang_mom = trexio.read_basis_shell_ang_mom(tf)
    shell_index = trexio.read_basis_shell_index(tf)
    exponent = trexio.read_basis_exponent(tf)
    coefficient = trexio.read_basis_coefficient(tf)
    nuc_charges = trexio.read_nucleus_charge(tf).astype(int)
    coords3d = trexio.read_nucleus_coord(tf)

    shells = list()
    for i in range(shell_num):
        mask = shell_index == i
        L = shell_ang_mom[i]
        nuc_index = nucleus_index[i]
        shell = Shell(
            L=L,
            atomic_num=nuc_charges[nuc_index],
            center=coords3d[nuc_index],
            coeffs=coefficient[mask],
            exps=exponent[mask],
            center_ind=nuc_index,
        )
        shells.append(shell)
    return TrexIOShells(shells, **kwargs)


def wavefunction_from_trexio(fn: str | Path, **kwargs):
    shell_kwargs = kwargs.pop("shell_kwargs", {})

    tf = trexio.File(str(fn), mode="r", back_end=trexio.TREXIO_HDF5)

    point_group = trexio.read_nucleus_point_group(tf)
    assert point_group == "C1"

    atoms = trexio.read_nucleus_label(tf)
    coords3d = trexio.read_nucleus_coord(tf)
    charges = trexio.read_nucleus_charge(tf)

    up_num = trexio.read_electron_up_num(tf)
    dn_num = trexio.read_electron_dn_num(tf)
    occ: tuple[int, int] = (up_num, dn_num)

    # Total spin
    Stot = 0.5 * up_num + -0.5 * dn_num
    mult = int(2 * Stot + 1)
    charge = sum(charges) - up_num - dn_num

    # ao_num, mo_num
    C_tio = trexio.read_mo_coefficient(tf).T
    mo_spin = trexio.read_mo_spin(tf)
    Ca = C_tio[:, mo_spin == 0]
    Cb = C_tio[:, mo_spin == 1]
    C = np.stack((Ca, Cb), axis=0)
    assert (len(Ca) + len(Cb)) == C_tio.shape[1]
    try:
        mo_ens = trexio.read_mo_energy(tf)
        mo_ens = mo_ens.reshape(2, -1)
    except trexio.Error:
        mo_ens = None

    restricted = up_num == dn_num
    # TODO: how to determine if calculation was done in an untrestricted formalism?
    # For now, check if the number of alpha and beta electrons differ and if
    # they are the same, check if the molecular orbtals differ.
    unrestricted = (not restricted) or (not np.allclose(Ca, Cb, atol=1e-10))

    bf_type = (
        BFType.CARTESIAN if trexio.read_ao_cartesian(tf) else BFType.PURE_SPHERICAL
    )
    ecp_electrons = None
    # Construct shells
    shells = shells_from_trexio(handle=tf, **shell_kwargs)
    tf.close()

    wf = Wavefunction(
        atoms=atoms,
        coords=coords3d.flatten(),
        charge=charge,
        mult=mult,
        unrestricted=unrestricted,
        occ=occ,
        C=C,
        bf_type=bf_type,
        shells=shells,
        ecp_electrons=ecp_electrons,
        mo_ens=mo_ens,
        **kwargs,
    )
    return wf


def wavefunction_to_trexio(
    wf: Wavefunction, fn: str | Path, point_group="C1", with_ovlp=False
):
    """Dump wavefunction to TREXIO file, overwriting an existing file."""
    fn = Path(fn)
    if fn.exists():
        fn.unlink()
    tf = trexio.File(str(fn), mode="w", back_end=trexio.TREXIO_HDF5)

    # 1. Metadata
    # Version seems to be written automatically
    # version = trexio.__version__
    # trexio.write_metadata_package_version(tf, version)
    trexio.write_metadata_code_num(tf, 1)
    trexio.write_metadata_code(tf, ["pysisyphus"])

    # 2. System
    #   2.1 Nucleus group
    nucleus_num = len(wf.atoms)
    trexio.write_nucleus_num(tf, nucleus_num)
    # TODO: ECP core electrons?!
    trexio.write_nucleus_charge(tf, wf.nuc_charges)
    trexio.write_nucleus_label(tf, wf.atoms)
    trexio.write_nucleus_point_group(tf, point_group)
    # trexio expects (3, nucleus.num) so we can use coords3d as it is
    trexio.write_nucleus_coord(tf, wf.coords3d)

    #   2.4 Electron group
    up_num, dn_num = wf.occ
    num = up_num + dn_num
    trexio.write_electron_num(tf, num)
    trexio.write_electron_up_num(tf, up_num)
    trexio.write_electron_dn_num(tf, dn_num)

    Ca, Cb = wf.C
    ao_num, mo_num = Ca.shape
    mo_spin = ([0] * mo_num) + ([1] * mo_num)
    # trexio expects MO coefficients (ao.num, mo.num) in column major order
    C = np.concatenate((Ca, Cb), axis=1).T
    mo_num = C.shape[0]
    trexio.write_ao_num(tf, ao_num)
    trexio.write_mo_num(tf, mo_num)
    # TODO: reorder
    is_cart = wf.is_cartesian
    # trexio uses the order 0, +1, -1, +2, -2, ..., +, -m
    # Matrix multiplying a pysisyphus' permutation matrix onto a vector holding
    # quantities in pysisyphus-order reorders to vector to the external order,
    # as defined in the permutation matrix.
    # /  -- pysis --> \
    # | |             |
    # | e             |
    # | x             |
    # | t             |
    # | ↓
    # \               /
    P = wf.get_permut_matrix()
    # TODO: also reorder energies
    C = np.einsum("ep,me->mp", P, C, optimize="greedy")
    if wf.mo_ens is not None:
        mo_ens = wf.mo_ens @ P
    else:
        mo_ens = None
    # If basis functions are Cartesian, then the alphabetical pysisyphus-ordering
    # is already correct.
    if not is_cart:
        trexio_shells = TrexIOShells.from_shells(wf.shells, ordering="native")
        P_trexio = trexio_shells.P_sph_native
        C = np.einsum("ep,mp->me", P_trexio, C, optimize="greedy")
        if mo_ens is not None:
            mo_ens = mo_ens @ P_trexio.T
    trexio.write_mo_coefficient(tf, C)
    if mo_ens is not None:
        trexio.write_mo_energy(tf, mo_ens.flatten())
    trexio.write_mo_spin(tf, mo_spin)

    # Write shells/basis functions/AOs
    cartesian = int(is_cart)
    trexio.write_ao_cartesian(tf, cartesian)

    # size nbfs
    ao_shells = list()
    ao_normalization = list()

    # size shell_num
    nucleus_index = list()
    shell_ang_mom = list()

    # size prim_num
    shell_index = list()
    exponent = list()
    coefficient = list()
    prim_factor = list()
    for i, shell in enumerate(wf.shells):
        size = shell.cartesian_size if cartesian else shell.sph_size
        ao_shells.extend([i] * size)

        nucleus_index.append(shell.center_ind)
        shell_ang_mom.append(shell.L)

        prim_size = len(shell.exps)
        shell_index.extend([i] * prim_size)
        exponent.extend(shell.exps.tolist())
        coefficient.extend(shell.coeffs_org.tolist())
        prim_factor.extend([norm_pgto(exp, shell.L) for exp in shell.exps])

    basis_type = "Gaussian"
    prim_num = len(exponent)
    shell_num = len(wf.shells)
    shell_factor = [1.0] * shell_num

    ao_normalization = np.ones_like(ao_shells)
    trexio.write_ao_shell(tf, ao_shells)
    trexio.write_ao_normalization(tf, ao_normalization)

    trexio.write_basis_type(tf, basis_type)
    trexio.write_basis_prim_num(tf, prim_num)
    trexio.write_basis_shell_num(tf, shell_num)
    trexio.write_basis_nucleus_index(tf, nucleus_index)
    trexio.write_basis_shell_ang_mom(tf, shell_ang_mom)
    trexio.write_basis_shell_factor(tf, shell_factor)
    trexio.write_basis_shell_index(tf, shell_index)
    trexio.write_basis_exponent(tf, exponent)
    trexio.write_basis_coefficient(tf, coefficient)
    trexio.write_basis_prim_factor(tf, prim_factor)

    if with_ovlp:
        S = wf.S
        if not is_cart:
            S = np.einsum("ep,ef,fq->pq", P, S, P, optimize="greedy")
            S = np.einsum("ep,pq,fq->ef", P_trexio, S, P_trexio, optimize="greedy")
            # TODO:: multiply P and P_trexio into one matrix
        trexio.write_ao_1e_int_overlap(tf, S)
    tf.close()
