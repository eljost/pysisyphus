#!/usr/bin/env python3

import argparse
import dataclasses
import enum
import functools
import hashlib
import operator
from pathlib import Path
import subprocess
import sys
from typing import Literal, Optional, Self, Sequence

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators import OverlapCalculator as ovlpcalc
from pysisyphus.calculators.ORCA import (
    parse_orca_cis,
    parse_orca_all_energies,
)
from pysisyphus.constants import AU2EV, BOHR2ANG
from pysisyphus.helpers import get_coords_diffs
import pysisyphus.io.cube as iocube
from pysisyphus.wavefunction import Shells, Wavefunction
from pysisyphus.wavefunction.diabatization import (
    dq_diabatization,
    edmiston_ruedenberg_diabatization_jacobi,
)
from pysisyphus.wavefunction.excited_states import (
    detachment_attachment_density,
    get_state_to_state_transition_density,
    make_density_matrices_for_root,
    norm_ci_coeffs,
)

# pysisyphus_addons imports below
from pysisyphus_addons.diabatization import intor
from pysisyphus_addons.grid import grid


class DiaKind(enum.Flag):
    NONE = 0
    EDMISTON_RUEDENBERG = enum.auto()
    BOYS = enum.auto()


DiaKindKeys = [key for key in DiaKind.__members__.keys() if key != "NONE"]


class CubeKind(enum.Flag):
    NONE = 0
    # Detachment & attachment densities
    ADIA_DA = enum.auto()
    DIA_DA = enum.auto()
    # Spin density
    ADIA_SD = enum.auto()
    DIA_SD = enum.auto()
    SD = ADIA_SD | DIA_SD


CubeKindKeys = [key for key in CubeKind.__members__.keys() if key != "NONE"]


def get_dia_result_fn(kind, base_name):
    return f"{base_name}_{kind}_dia_result.npz"


def get_states_hash(states: np.ndarray):
    inp_str = "-".join(map(str, states))
    hash_object = hashlib.md5(inp_str.encode("utf-8"))
    return hash_object.hexdigest()


@dataclasses.dataclass
class DiaInput:
    wf: Wavefunction
    all_ens: np.ndarray
    Xa: np.ndarray
    Ya: np.ndarray
    Xb: np.ndarray
    Yb: np.ndarray
    states: np.ndarray
    base_name: str

    @property
    def adia_ens(self):
        ae = self.all_ens
        # Set GS to 0.0
        return (ae - ae[0])[self.states]

    @property
    def Ca(self):
        return self.wf.C[0]

    @property
    def Cb(self):
        return self.wf.C[1]

    def update_states(self, other: Self | None):
        new_states = get_new_states(self, other)
        self.states = new_states


def get_new_states(
    dia_inp: DiaInput, dia_inp_prev: Optional[DiaInput] = None, ovlp_type="tden"
):
    # TODO: add logging
    if dia_inp_prev is None:
        return dia_inp.states

    S_AO = ovlpcalc.get_sao_from_mo_coeffs(dia_inp.Ca)
    di = dia_inp
    dip = dia_inp_prev
    ovlp_mat = ovlpcalc.get_ovlp_mat(
        dip.Ca,
        dip.Cb,
        dip.Xa,
        dip.Ya,
        dip.Xb,
        dip.Yb,
        di.Ca,
        di.Cb,
        di.Xa,
        di.Ya,
        di.Xb,
        di.Yb,
        S_AO,
        ovlp_type,
    )
    states = np.zeros_like(dip.states)
    for i, root_prev in enumerate(dip.states):
        if root_prev == 0:
            root_new = 0
            msg = f"Keeping root {root_prev}"
        else:
            root_new = ovlpcalc.root2_from_ovlp_mat_and_root1(ovlp_mat, root_prev)
            ovlp = ovlp_mat[root_prev - 1, root_new - 1]
            msg = (
                f"Overlap between previous root {root_prev} "
                f"and new root {root_new} is {ovlp: >6.2%}"
            )
        print(msg)
        states[i] = root_new
    return states


def row_major_iter(nitems: int):
    """Row-major iterator over pairs indices."""
    for i in range(nitems):
        for j in range(i + 1):
            yield i, j


def map2d(i, j):
    """Map 2d index to 1d w/ row-major ordering."""
    if j > i:
        i, j = j, i
    return i * (i + 1) // 2 + j


def make_ao_spin_densities(
    C: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    states: Sequence[int],
    spin: Literal["alpha", "beta"],
    P: Optional[np.ndarray] = None,
):
    """Construct alpha or beta spin (transition) densities in the AO basis.

    Paramters
    ---------
    C
        Either alpha or beta MO coefficient matrix. 2d array of shape (naos, nmos).
    X
        X vector as obtained from a TD-DFT/TDA calculation. 3d array of shape
        (nstates, nocc, nvirt).
    Y
        Y vector as obtained from a TD-DFT/TDA calculation. 3d array of shape
        (nstates, nocc, nvirt). Will be a zero array in TDA/CIS calculations.
    states
        Sequence of integers containing state indices for which the (transition)
        densities will be formed. Indices in this sequence should be consisent
        with the densities in X and Y. If one wants to stick with actual state
        indices (0, 1, 2, ...) then X and Y should contain (dummy) ground state
        transition densities that are zero throughout.
    spin
        Either "alpha" or "beta".
    P
        Optional permutation matrix that is used to bring the basi functions (AOs)
        to pysisyphus-ordering. If given, P is a 2d permutation matrix of shape
        (naos, naos) with the external order along the rows and pysisyphus order
        along the columns.

    Returns
    -------
    densities
        3d array of shape (nstates * (nstates + 1) // 2, naos, naos) containing
        AO (transition) densities with AOs in pysisyphus order.
    state_pairs
        List of integer tuples containing the state pair indices.
    """
    nstates = len(states)

    _, nocc, nvirt = X.shape
    naos, nmos = C.shape
    # Quick sanity check if X and C are consistent
    assert nocc + nvirt == nmos

    Cocc = C[:, :nocc]
    Cvirt = C[:, nocc:]
    # Occ-occ part of the ground state in the MO basis.
    eye = np.eye(nocc)

    ndens = nstates * (nstates + 1) // 2
    densities = np.zeros((ndens, naos, naos))
    state_pairs = list()
    k = 0
    for i, I in enumerate(states):
        XI = X[I]
        YI = Y[I]
        for J in states[: i + 1]:
            assert I >= J
            state_pairs.append((I, J))
            # State density (I == J) or excited state to excited state density
            if I == J or J > 0:
                dens_mo = get_state_to_state_transition_density(
                    X[J], XI
                ) + get_state_to_state_transition_density(Y[J], YI)
                if I == J:
                    dens_mo[:nocc, :nocc] += eye
                dens_ao = C @ dens_mo @ C.T
            # Ground state to excited state. In this case we use the transition
            # densities directly.
            else:
                dens_ao = Cocc @ (XI + YI) @ Cvirt.T
            densities[k] = dens_ao
            print(f"Set ({I: >3d},{J: >3d}) {spin: >5s} spin density")
            k += 1

    # Convert to pysisyphus-order if permutation matrix was given
    if P is not None:
        densities = np.einsum("rm,sn,Irs->Imn", P, P, densities, optimize="greedy")
    return densities, state_pairs


def run_er_dia(
    wf: Wavefunction,
    densities: np.ndarray,
    adia_ens: np.ndarray,
    tensor_fn: Optional[Path] = None,
    aux_shells: Optional[Shells] = None,
):
    """Edmiston-Ruedenberg diabatization wrapper.

    Parameters
    ----------
    wf
        Pysisyphus wavefunction.
    densities
        3d array of shape (nstates * (nstates + 1) // 2, naos, naos)
        (transition) densities in the AO basis in pysisyphus order.
    adia_ens
        Adiabatic excitation energies in Hartee.
    aux_shells
        Shells object holding the auxiliary basis. If not provided, pysisyphus
        will use the universal Weigend def2-J aux basis.

    Returns
    -------
    dia_result
        DiabatizationResult object representing the results of the ER diabatization.
    """
    shells = wf.shells
    if aux_shells is None:
        # TODO: handle storage of aux basis
        aux_basis_fn = "def2-universal-jfit.json"
        aux_shells = shells.from_basis(aux_basis_fn)
    adia_ens_eV = adia_ens * AU2EV
    if tensor_fn is not None and tensor_fn.exists():
        print(f"Found existing Coulomb tensor in {tensor_fn}.")
        coulomb_tensor = np.load(tensor_fn)
    else:
        # This is were the magic happens and we call into Fortran
        coulomb_tensor = intor.contract_coulomb_densities_4d(
            shells, aux_shells, densities
        )
        # Skip saving when no fn was given
        if tensor_fn is not None:
            np.save(tensor_fn, coulomb_tensor)
            print(f"Saved Coulomb tensor to {tensor_fn}.")
    dia_result = edmiston_ruedenberg_diabatization_jacobi(adia_ens_eV, coulomb_tensor)
    return dia_result


def fmt_dpm(dpm: np.ndarray) -> str:
    """Format dipole moment."""
    return np.array2string(dpm, formatter={"float": lambda flt: f"{flt: >10.6f}"})


def make_dip_moms_2d(wf: Wavefunction, densities: np.ndarray):
    assert densities.ndim == 3
    kind = "coc"
    origin = wf.get_origin(kind)
    origin_ang = origin * BOHR2ANG
    print(f"Using '{kind}'-origin at {origin} au ({origin_ang} Å).")
    # Calcualte dipole moment integrals in AO basis
    dip_ints = wf.dipole_ints(origin)
    print("Calculated dipole integrals")
    # Bring dipole moment integrals into pysisyphus order
    P = wf.get_permut_matrix()
    dip_ints = np.einsum("Irs,rm,sn->Imn", dip_ints, P, P, optimize="greedy")
    dip_moms_2d = np.einsum("imn,Kmn->Ki", dip_ints, densities, optimize="greedy")
    return dip_moms_2d


def run_boys_dia(dip_moms_2d: np.ndarray, adia_ens: np.ndarray):
    """Boys-diabatization wrapper.

    Parameters
    ----------
    dip_moms_2d
        2d array of shape (nstates * (nstates + 1) // 2, 3) containing
        (transition) dipole moments in row-major order.
    adia_ens
        Adiabatic excitation energies in Hartee.

    Returns
    -------
    dia_result
        DiabatizationResult object representing the results of the Boys-diabatization.
    """

    adia_ens_eV = adia_ens * AU2EV
    nstates = len(adia_ens)
    dip_moms = np.zeros((3, nstates, nstates))
    # Bring dipole moments into a shape suitable for the diabatization procedure
    for (i, j), dpm in zip(row_major_iter(nstates), dip_moms_2d, strict=True):
        dip_moms[:, i, j] = dpm
        dip_moms[:, j, i] = dpm
    dia_result = dq_diabatization(adia_ens_eV, dip_moms)
    return dia_result


def make_ao_da_densities(
    wf: Wavefunction,
    roots: Sequence[int],
    Xa: np.ndarray,
    Ya: np.ndarray,
    Xb: np.ndarray,
    Yb: np.ndarray,
):
    """Calculate detachment/attachment densities in the AO basis."""

    # Any root <= 0 will yield erronous results when rootm1 is calculated, as
    # root 0 would become -1, which would refer to the last state!
    assert min(roots) > 0, "ES roots must be > 0!"
    max_root = max(roots)
    nstates = len(Xa)
    assert max_root <= nstates

    nroots = len(roots)
    # In AO basis
    Pa, Pb = wf.P
    # Overlap matrix
    S = wf.S
    # Transform density matrix to MO basis
    Ca, Cb = wf.C
    Pa_mo = Ca.T @ S @ Pa @ S @ Ca
    Pb_mo = Cb.T @ S @ Pb @ S @ Cb
    print("Formed GS density matrix in the MO basis")

    det_dens_mo_a = np.zeros((nroots, *Pa_mo.shape))
    att_dens_mo_a = np.zeros_like(det_dens_mo_a)
    det_dens_mo_b = np.zeros((nroots, *Pb_mo.shape))
    att_dens_mo_b = np.zeros_like(det_dens_mo_b)
    # TODO: take care of GS that has norm 0.0
    for i, root in enumerate(roots):
        rootm1 = root - 1
        # In MO basis
        Pa_exc, Pb_exc = make_density_matrices_for_root(
            rootm1,
            wf.restricted,
            Xa,
            Ya,
            Xb,
            Yb,
        )

        diff_dens_a = Pa_exc - Pa_mo
        P_det_a, P_att_a = detachment_attachment_density(diff_dens_a, verbose=True)
        det_dens_mo_a[i] = P_det_a
        att_dens_mo_a[i] = P_att_a

        diff_dens_b = Pb_exc - Pb_mo
        P_det_b, P_att_b = detachment_attachment_density(diff_dens_b, verbose=True)
        det_dens_mo_b[i] = P_det_b
        att_dens_mo_b[i] = P_att_b

        print(f"Calculated A&D density matrices in the MO basis for {root=}")

    def transform_mo_ao(C, dens_mo):
        return np.einsum("mr,Irs,ns->Imn", C, dens_mo, C, optimize="greedy")

    det_dens_ao_a = transform_mo_ao(Ca, det_dens_mo_a)
    att_dens_ao_a = transform_mo_ao(Ca, att_dens_mo_a)
    det_dens_ao_b = transform_mo_ao(Cb, det_dens_mo_b)
    att_dens_ao_b = transform_mo_ao(Cb, att_dens_mo_b)
    return det_dens_ao_a, att_dens_ao_a, det_dens_ao_b, att_dens_ao_b


def cube_densities(
    wf, densities, labels, out_dir: Path, reorder_coeffs=None, **grid_kwargs
):
    _grid_kwargs = {
        "num": 50,
        "margin": 4.0,
    }
    _grid_kwargs.update(grid_kwargs)

    # When called w/o reorder_coeffs, we expect densities over spherical
    # basis functions in external-order. Then we transform them to densities
    # over Cartesian basis functions in pysisyphus-order.
    if reorder_coeffs is None:
        reorder_coeffs = wf.shells.reorder_c2s_coeffs
    densities_cart = np.einsum(
        "Jmn,mo,np->Jop", densities, reorder_coeffs, reorder_coeffs, optimize="greedy"
    )

    ndens = len(densities_cart)
    grid3d, spacing, grid_shape = iocube.get_grid(wf, **_grid_kwargs)
    axes = np.diag(spacing)
    origin = grid3d[0]
    vol_elem = np.prod(spacing)

    # Calculate densities on the grid using pysisyphus-addons
    densities_grid = grid.eval_densities(wf.shells, grid3d, densities_cart)
    vol_datas = densities_grid.reshape(ndens, *grid_shape)
    cube_fns = list()
    vol_fmt = " >8.3e"
    for vol_data, label in zip(vol_datas, labels):
        # Calculate number of electrons from the density.
        N = vol_data.sum() * vol_elem
        cube = iocube.Cube.from_wf_and_grid(
            wf, vol_data.reshape(grid_shape), origin, axes
        )
        cube_fn = out_dir / f"{label}.cub"
        cube.write(cube_fn)
        cube_fns.append(cube_fn)
        print(
            f"{cube_fn}: min={vol_data.min():{vol_fmt}}, max={vol_data.max():{vol_fmt}}, "
            f"{N=: >8.4f}, {grid_shape=}"
        )
    return cube_fns


def cube_da_densities(
    wf, states, Xa, Ya, Xb, Yb, kind: str, base: str, out_dir: Path, **grid_kwargs
):
    nstates = len(states)
    det_a, att_a, det_b, att_b = make_ao_da_densities(wf, states, Xa, Ya, Xb, Yb)
    # Total detachment and attachment densities
    det_dens = det_a + det_b
    att_dens = att_a + att_b
    da_densities = np.concatenate((det_dens, att_dens), axis=0)
    # Construct cube labels
    kinds = [f"{kind}_det"] * nstates + [f"{kind}_att"] * nstates
    label_states = list(states)
    labels = [
        f"{base}_{kind}{state}"
        for state, kind in zip(label_states + label_states, kinds)
    ]
    cube_fns = cube_densities(wf, da_densities, labels, out_dir, **grid_kwargs)
    return cube_fns


def make_diabatic_properties(U, adia_props: np.ndarray):
    assert adia_props.ndim > 2
    nadia, ndia = U.shape
    assert nadia == ndia
    assert adia_props.shape[0] == adia_props.shape[1] == nadia

    add_shape = [1] * (adia_props.ndim - 2)
    dia_props = np.zeros_like(adia_props)
    for a, A in enumerate(U.T):
        for b, B in enumerate(U.T[a:], a):
            mix = np.outer(A, B)
            mix = mix.reshape(*mix.shape, *add_shape)
            dia_props[a, b] = (mix * adia_props).sum(axis=(0, 1))
            if a != b:
                dia_props[b, a] = dia_props[a, b]
    return dia_props


@functools.singledispatch
def run_dia(
    wf: Wavefunction,
    adia_ens: np.ndarray,
    Xa: np.ndarray,
    Ya: np.ndarray,
    Xb: np.ndarray,
    Yb: np.ndarray,
    states: np.ndarray | Sequence[int],
    base_name: str,
    dia_kinds: DiaKind,
    cube_kinds: CubeKind,
    grid_kwargs: Optional[dict] = None,
    out_dir: Optional[Path] = None,
    force: bool = False,
):
    """
    Parameters
    ----------
    force
        Force calculation of the Coulomb tensor and don't try to save it to or
        load it from disk. Useful for testing.
    """
    if grid_kwargs is None:
        grid_kwargs = {}

    if out_dir is None:
        out_dir = Path(".")
    if not out_dir.exists():
        out_dir.mkdir()

    cube_da_densities_wrapped = functools.partial(
        cube_da_densities, out_dir=out_dir, base_name=base_name, **grid_kwargs
    )

    states = np.sort(states, kind="stable")
    nstates = len(states)

    Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)

    # Add fake/empty zero ground state transition density
    # Alpha
    _, nocca, nvirta = Xa.shape
    gs_zeros_a = np.zeros((1, nocca, nvirta))
    Xa = np.concatenate((gs_zeros_a, Xa), axis=0)
    Ya = np.concatenate((gs_zeros_a, Ya), axis=0)
    # Beta
    _, noccb, nvirtb = Xb.shape
    gs_zeros_b = np.zeros((1, noccb, nvirtb))
    Xb = np.concatenate((gs_zeros_b, Xb), axis=0)
    Yb = np.concatenate((gs_zeros_b, Yb), axis=0)

    # Construct (transition) densities for the states in the diabatization.
    # The densities will be used to evaluated properites (Coulomb tensor, multipole moments).
    Ca, Cb = wf.C
    naos, _ = Ca.shape

    P = wf.get_permut_matrix()
    adia_densities_a, state_pairs = make_ao_spin_densities(
        Ca, Xa, Ya, states, "alpha", P=P
    )
    adia_densities_b, _ = make_ao_spin_densities(Cb, Xb, Yb, states, "beta", P=P)
    adia_densities = adia_densities_a + adia_densities_b

    nstates = len(adia_ens)
    # Report adiabatic (transiton) dipole moments calculated with our densities.
    # These quantities are often calculated by the underlying QC codes themselves,
    # so the values printed here can be compared to the QC code output.
    #
    # The dipole integrals will also be later passed to the Boys diabatization routine,
    # if enabled.
    dip_moms_2d = make_dip_moms_2d(wf, adia_densities)
    print("(Transition) dipole moments:")
    # TODO: also report GS-> ES DPMs as these are usually printed by the QC codes
    for (i, j), dpm in zip(row_major_iter(nstates), dip_moms_2d, strict=True):
        I = states[i]
        J = states[j]
        print(f"\t{I: >03d} -> {J: >03d}: {fmt_dpm(dpm)} au")

    cube_fns = {}
    dia_results = {}

    # Adiabatic detachment/attachment densities
    if CubeKind.ADIA_DA in cube_kinds:
        # Skip GS in adiabatic detachment/attachment calculation.
        states_no_gs = states[states > 0]
        adia_da_cube_fns = cube_da_densities_wrapped(
            wf,
            states_no_gs,
            # Don't pass empty GS arrays
            Xa[1:],
            Ya[1:],
            Xb[1:],
            Yb[1:],
            "adia",
            # base_name,
            # out_dir,
            # **grid_kwargs,
        )
        cube_fns["ADIA_DA"] = adia_da_cube_fns

    # Calculate spin densities if spin density cubes were requested
    if cube_kinds & CubeKind.SD:
        c2s_coeffs = wf.shells.cart2sph_coeffs
        adia_spin_densities_3d = adia_densities_a - adia_densities_b
        adia_spin_densities = np.zeros((nstates, nstates, naos, naos))
        for (I, J), sd in zip(
            row_major_iter(nstates), adia_spin_densities_3d, strict=True
        ):
            adia_spin_densities[I, J] = sd
            if I != J:
                adia_spin_densities[J, I] = sd

    # Adiabatic spin densities
    if CubeKind.ADIA_SD in cube_kinds:
        adia_state_spin_densities = np.array(
            [adia_spin_densities[i, i] for i in range(nstates)]
        )
        adia_state_sd_labels = [f"{base_name}_adia_sd{I}" for I in states]
        cube_densities(
            wf,
            adia_state_spin_densities,
            adia_state_sd_labels,
            reorder_coeffs=c2s_coeffs,
            **grid_kwargs,
        )

    #########################
    # Actual diabatizations #
    #########################

    add_savez_kwargs = {
        "states": states,
    }

    get_result_fn = functools.partial(get_dia_result_fn, base_name=base_name)

    # Edmiston-Ruedenberg diabatization
    if DiaKind.EDMISTON_RUEDENBERG in dia_kinds:
        states_hash = get_states_hash(states)
        if force:
            tensor_fn = None
        else:
            tensor_fn = out_dir / Path(f"{base_name}_coulomb_tensor_{states_hash}.npy")
        er_dia_result = run_er_dia(wf, adia_densities, adia_ens, tensor_fn)
        er_dia_results_fn = get_result_fn("er")
        er_dia_result.savez(out_dir / er_dia_results_fn, **add_savez_kwargs)
        dia_results["er"] = er_dia_result
        print(er_dia_result.render_report())
        sys.stdout.flush()

    # Boys diabatization
    if DiaKind.BOYS in dia_kinds:
        boys_dia_result = run_boys_dia(dip_moms_2d, adia_ens)
        boys_dia_results_fn = get_result_fn("boys")
        boys_dia_result.savez(out_dir / boys_dia_results_fn, **add_savez_kwargs)
        dia_results["boys"] = boys_dia_result
        print(boys_dia_result.render_report())
        sys.stdout.flush()

    if CubeKind.DIA_DA in cube_kinds and 0 in states:
        cube_kinds &= ~CubeKind.DIA_DA
        print(
            "Ground state involved! Disabled generation of diabatic "
            "detachment/attachment densities!"
        )

    # Loop over all diabatization that were carried out and calculate diabatic
    # properties as requested.
    for prefix, dia_result in dia_results.items():
        # One column per diabatic state
        U = dia_result.U
        Xa_dia = Xa[states].copy()
        Xa_dia = np.einsum("JI,Jia->Iia", U, Xa_dia)
        Ya_dia = Ya[states].copy()
        Ya_dia = np.einsum("JI,Jia->Iia", U, Ya_dia)
        Xb_dia = Xb[states].copy()
        Xb_dia = np.einsum("JI,Jia->Iia", U, Xb_dia)
        Yb_dia = Yb[states].copy()
        Yb_dia = np.einsum("JI,Jia->Iia", U, Yb_dia)
        # Renumber the adiabatic states
        dia_states = (np.arange(nstates) + 1).tolist()

        # Diabatic detachment/attachment densities
        if CubeKind.DIA_DA in cube_kinds:
            assert 0 not in states
            dia_da_cube_fns = cube_da_densities_wrapped(
                wf,
                dia_states,
                Xa_dia,
                Ya_dia,
                Xb_dia,
                Yb_dia,
                f"{prefix}_dia",
                # base_name,
                # out_dir,
                # **grid_kwargs,
            )
            cube_fns[f"{prefix.upper()}_DIA_DA"] = dia_da_cube_fns

        # Diabatic spin densities
        if CubeKind.DIA_SD in cube_kinds:
            dia_spin_densities = make_diabatic_properties(U, adia_spin_densities)
            dia_state_spin_densities = np.array(
                [dia_spin_densities[i, i] for i in range(nstates)]
            )
            dia_state_sd_labels = [
                f"{base_name}_{prefix}_dia_sd{I}" for I in dia_states
            ]
            dia_sd_cube_fns = cube_densities(
                wf,
                dia_state_spin_densities,
                dia_state_sd_labels,
                reorder_coeffs=c2s_coeffs,
                **grid_kwargs,
            )
            cube_fns[f"{prefix.upper()}_DIA_SD"] = dia_sd_cube_fns

    return dia_results, cube_fns


@run_dia.register
def _(dia_inp: DiaInput, **kwargs):
    return run_dia(
        dia_inp.wf,
        adia_ens=dia_inp.adia_ens,
        Xa=dia_inp.Xa,
        Ya=dia_inp.Ya,
        Xb=dia_inp.Xb,
        Yb=dia_inp.Yb,
        states=dia_inp.states,
        base_name=dia_inp.base_name,
        **kwargs,
    )


def parse_orca(base: str | Path, triplets: bool = False):
    base = Path(base)
    wf_fn = base.with_suffix(".bson")
    if not wf_fn.exists():
        subprocess.run(f"orca_2json {str(base)} -bson", shell=True)
    wf = Wavefunction.from_file(wf_fn)
    print(f"Read wavefunction: {wf}")

    cis_fn = base.with_suffix(".cis")
    log_fn = base.with_suffix(".out")
    if not log_fn.exists():
        log_fn = base.with_suffix(".log")

    all_ens = parse_orca_all_energies(log_fn, do_tddft=True, triplets=triplets)
    print("Parsed energies from ORCA log file")

    Xa, Ya, Xb, Yb = parse_orca_cis(
        cis_fn, restricted_same_ab=True, triplets_only=triplets
    )
    print("Parsed ORCA .cis file")
    return wf, all_ens, Xa, Ya, Xb, Yb


def plot_dia_results(coord_diffs, dia_results, states, kind):
    assert kind in ("boys", "er")
    all_adia_ens = list()
    all_couplings = list()
    for dia_res in dia_results:
        adia_ens = dia_res.adia_ens
        all_adia_ens.append(adia_ens)
        coupling01 = dia_res.dia_mat[0, 1]
        # TODO: way to get unique couplings
        all_couplings.append(coupling01)

    all_adia_ens = np.array(all_adia_ens)
    all_couplings = np.array(all_couplings)
    abs_all_couplings = np.abs(all_couplings)

    fig, ax = plt.subplots()
    for state, state_ens in zip(states, all_adia_ens.T):
        ax.plot(coord_diffs, state_ens, "o-", label=state)
    ax.set_xlabel("Interpolation coordinate")
    ax.set_ylabel("Adiabatic energies / eV")
    ax.legend()

    ax2 = ax.twinx()
    ax2.plot(coord_diffs, abs_all_couplings, "o-", c="green", label="Coupling")
    ax2.set_ylabel("abs. diabatic coupling / eV")
    ax.set_title(f"{kind.capitalize()}-Diabatization")

    ax.set_xlim(coord_diffs[0], coord_diffs[-1])
    fig.tight_layout()
    fig.savefig(f"{kind}_diabatization.pdf")
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("orca_outs", nargs="+")
    parser.add_argument("--states", nargs="+", type=int, required=True)
    parser.add_argument("--triplets", action="store_true")
    parser.add_argument(
        "--ovlp",
        action="store_true",
        help="Determine states from tden-overlaps between steps.",
    )
    # Diabatization algorithm selection
    parser.add_argument(
        "--dia",
        type=lambda s: DiaKind[s.upper()],
        default=[DiaKind.NONE],
        nargs="+",
        help=f"Available diabatization algorithms: {', '.join(DiaKindKeys)}.",
    )
    # Cube generation selection
    parser.add_argument(
        "--cube",
        type=lambda s: CubeKind[s.upper()],
        default=[CubeKind.NONE],
        nargs="+",
        help=(
            "Available cubes (DA: detachment/attachment, SD: spin density): "
            f"{', '.join(CubeKindKeys)}."
        ),
    )
    parser.add_argument(
        "--grid-points",
        type=int,
        default=50,
        help="Number of grid points per cube axis.",
    )

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])

    orca_outs = list(map(Path, args.orca_outs))
    states = args.states
    triplets = args.triplets
    dia_kinds = functools.reduce(operator.or_, args.dia)
    cube_kinds = functools.reduce(operator.or_, args.cube)
    grid_points = args.grid_points
    ovlp = args.ovlp

    grid_kwargs = {
        "num": grid_points,
    }

    nouts = len(orca_outs)
    nstates = len(states)

    all_dia_results = list()
    all_coords = list()
    dia_inp_prev = None
    all_states = np.zeros((nouts, nstates))
    for i, base in enumerate(orca_outs):
        base_name = base.stem
        wf, all_ens, Xa, Ya, Xb, Yb = parse_orca(base, triplets=triplets)
        sys.stdout.flush()

        dia_inp = DiaInput(
            wf=wf,
            all_ens=all_ens,
            Xa=Xa,
            Ya=Ya,
            Xb=Xb,
            Yb=Yb,
            states=states,
            base_name=base_name,
        )
        # Update states by calculating overlaps if requested
        if ovlp:
            # TODO: report previous & new state before/after overlaps
            dia_inp.update_states(dia_inp_prev)
        print(f"States for {base} are {dia_inp.states}")
        all_states[i] = states

        dia_results, cube_fns = run_dia(
            dia_inp,
            dia_kinds=dia_kinds,
            cube_kinds=cube_kinds,
            grid_kwargs=grid_kwargs,
        )
        all_dia_results.append(dia_results)
        all_coords.append(wf.coords)
        dia_inp_prev = dia_inp

    np.save("all_states.npy", all_states)

    if nouts > 1:
        cds = get_coords_diffs(all_coords, align=True)
        if DiaKind.EDMISTON_RUEDENBERG in dia_kinds:
            er_dia_results = [dr["er"] for dr in all_dia_results]
            plot_dia_results(cds, er_dia_results, states, "er")
        if DiaKind.BOYS in dia_kinds:
            boys_dia_results = [dr["boys"] for dr in all_dia_results]
            plot_dia_results(cds, boys_dia_results, states, "boys")


if __name__ == "__main__":
    main()
