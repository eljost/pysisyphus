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
import pysisyphus.io.hdf5 as ioh5
from pysisyphus.wavefunction import Shells, Wavefunction
from pysisyphus.diabatization import logger
from pysisyphus.diabatization import plot as dia_plot, chimerabatch
from pysisyphus.diabatization.coulomb import edmiston_ruedenberg_diabatization_jacobi
from pysisyphus.diabatization.coulomb_eta import edmiston_ruedenberg_eta_diabatization
from pysisyphus.diabatization.multipole import dq_diabatization
from pysisyphus.wavefunction.excited_states import (
    detachment_attachment_density,
    get_state_to_state_transition_density,
    make_density_matrices_for_root,
    norm_ci_coeffs,
)


try:
    # pysisyphus_addons imports below
    from pysisyphus_addons.diabatization import intor
    from pysisyphus_addons.grid import grid

    can_pysis_addons = True
except (ImportError, ModuleNotFoundError):
    can_pysis_addons = False


class DiaKind(enum.Flag):
    NONE = 0
    EDMISTON_RUEDENBERG = enum.auto()
    EDMISTON_RUEDENBERG_ETA = enum.auto()
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


def get_dia_data_model(kind, nstates):
    data_model = {
        "states": (nstates,),
        "U": (nstates, nstates),
        "adia_ens": (nstates,),
    }
    if kind in ("er", "ereta"):
        data_model.update(
            {
                "R_tensor": (nstates, nstates, nstates, nstates),
            }
        )
    if kind in ("boys",):
        data_model.update(
            {
                "dip_mom_tensor": (3, nstates, nstates),
            }
        )
    return data_model


def get_states_hash(states: np.ndarray):
    inp_str = "-".join(map(str, states))
    hash_object = hashlib.md5(inp_str.encode("utf-8"))
    return hash_object.hexdigest()


def plot_diabatic(dia_res, key, states, base_name, out_dir):
    G = dia_plot.state_graph_from_en_mat(dia_res.dia_mat, state_inds=states)
    fig = dia_plot.draw_state_graph(G)
    fig_fn = out_dir / f"{base_name}_{key}_plot.pdf"
    fig.savefig(fig_fn)
    logger.info(f"Wrote plot of diabatic states to '{fig_fn}'.")


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
        logger.info(msg)
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
    full: bool = True,
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
    full
        Whether to include density terms arising from occupied orbitals in the
        ground state. For Edmiston-Ruedenberg-Eta diabatization this term has
        to be excluded.

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
                if full and I == J:
                    dens_mo[:nocc, :nocc] += eye
                dens_ao = C @ dens_mo @ C.T
            # Ground state to excited state. In this case we use the transition
            # densities directly.
            else:
                dens_ao = Cocc @ (XI + YI) @ Cvirt.T
            densities[k] = dens_ao
            logger.info(f"Set ({I: >3d},{J: >3d}) {spin: >5s} spin density")
            k += 1

    # Convert to pysisyphus-order if permutation matrix was given
    if P is not None:
        densities = np.einsum("rm,sn,Irs->Imn", P, P, densities, optimize="greedy")
    return densities, state_pairs


def run_coulomb_dia(
    wf: Wavefunction,
    densities: np.ndarray,
    adia_ens: np.ndarray,
    dia_func,
    tensor_fn: Optional[Path] = None,
    aux_shells: Optional[Shells] = None,
    **dia_kwargs,
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
        logger.info(f"Found existing Coulomb tensor in {tensor_fn}.")
        coulomb_tensor = np.load(tensor_fn)
    else:
        # This is were the magic happens and we call into Fortran
        coulomb_tensor = intor.contract_coulomb_densities_4d(
            shells, aux_shells, densities
        )
        # Skip saving when no fn was given
        if tensor_fn is not None:
            np.save(tensor_fn, coulomb_tensor)
            logger.info(f"Saved Coulomb tensor to {tensor_fn}.")
    dia_result = dia_func(adia_ens_eV, coulomb_tensor, **dia_kwargs)
    return dia_result


def fmt_dpm(dpm: np.ndarray) -> str:
    """Format dipole moment."""
    return np.array2string(dpm, formatter={"float": lambda flt: f"{flt: >10.6f}"})


def make_dip_moms_2d(wf: Wavefunction, densities: np.ndarray):
    assert densities.ndim == 3
    kind = "coc"
    origin = wf.get_origin(kind)
    origin_ang = origin * BOHR2ANG
    logger.info(f"Using '{kind}'-origin at {origin} au ({origin_ang} Å).")
    # Calcualte dipole moment integrals in AO basis
    dip_ints = wf.dipole_ints(origin)
    logger.info("Calculated dipole integrals")
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
    logger.info("Formed GS density matrix in the MO basis")

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

        logger.info(f"Calculated A&D density matrices in the MO basis for {root=}")

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
        logger.info(
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


def postprocess_dia(
    kind,
    dia_result,
    out_dir,
    h5_fn,
    plot_diabatic,
    all_dia_results,
    add_attrs: Optional[dict] = None,
    add_datasets: Optional[dict] = None,
):
    if add_attrs is None:
        add_attrs = {}

    h5_fn = out_dir / h5_fn
    h5_data_model = get_dia_data_model(kind, dia_result.nstates)
    h5_group = ioh5.get_h5_group(h5_fn, kind, data_model=h5_data_model, reset=True)
    h5_group.attrs["kind"] = kind
    for h5_key in h5_data_model:
        try:
            val = getattr(dia_result, h5_key)
        except AttributeError:
            val = add_datasets[h5_key]
        h5_group[h5_key][:] = val
    for k, v in add_attrs.items():
        h5_group[k] = v

    all_dia_results[kind] = dia_result
    logger.info(dia_result.render_report())
    plot_diabatic(dia_result, kind)
    sys.stdout.flush()


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
    dia_kwargs: Optional[dict] = None,
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
    if dia_kwargs is None:
        dia_kwargs = {}
    if grid_kwargs is None:
        grid_kwargs = {}

    if out_dir is None:
        out_dir = Path(".")
    if not out_dir.exists():
        out_dir.mkdir()

    all_dia_results = {}

    cube_da_densities_wrapped = functools.partial(
        # cube_da_densities, out_dir=out_dir, base_name=base_name, **grid_kwargs
        cube_da_densities,
        out_dir=out_dir,
        base=base_name,
        **grid_kwargs,
    )
    plot_diabatic_wrapped = functools.partial(
        plot_diabatic, states=states, base_name=base_name, out_dir=out_dir
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
    # See eq. (A5) in [???] and the comment below. For Edmiston-Ruedenberg-Eta
    # diabatization we have to exclude the diagonal terms.
    full = DiaKind.EDMISTON_RUEDENBERG_ETA not in dia_kinds
    adia_densities_a, _ = make_ao_spin_densities(
        Ca, Xa, Ya, states, "alpha", P=P, full=full
    )
    adia_densities_b, _ = make_ao_spin_densities(
        Cb, Xb, Yb, states, "beta", P=P, full=full
    )
    adia_densities = adia_densities_a + adia_densities_b

    nstates = len(adia_ens)
    # Report adiabatic (transiton) dipole moments calculated with our densities.
    # These quantities are often calculated by the underlying QC codes themselves,
    # so the values printed here can be compared to the QC code output.
    #
    # The dipole integrals will also be later passed to the Boys diabatization routine,
    # if enabled.
    dip_moms_2d = make_dip_moms_2d(wf, adia_densities)
    logger.info("(Transition) dipole moments:")
    # TODO: also report GS-> ES DPMs as these are usually printed by the QC codes
    for (i, j), dpm in zip(row_major_iter(nstates), dip_moms_2d, strict=True):
        I = states[i]
        J = states[j]
        logger.info(f"\t{I: >03d} -> {J: >03d}: {fmt_dpm(dpm)} au")

    add_datasets = {
        "states": states,
    }
    h5_fn = f"{base_name}_dia_result.h5"
    postprocess_dia_wrapped = functools.partial(
        postprocess_dia,
        out_dir=out_dir,
        h5_fn=h5_fn,
        plot_diabatic=plot_diabatic_wrapped,
        all_dia_results=all_dia_results,
        add_datasets=add_datasets,
    )

    cube_fns = {}
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
            out_dir,
            reorder_coeffs=c2s_coeffs,
            **grid_kwargs,
        )

    #########################
    # Actual diabatizations #
    #########################

    states_hash = get_states_hash(states)
    if force:
        coulomb_tensor_fn = None
    else:
        coulomb_tensor_fn = out_dir / Path(
            f"{base_name}_coulomb_tensor_{states_hash}.npy"
        )

    # Edmiston-Ruedenberg diabatization
    if DiaKind.EDMISTON_RUEDENBERG in dia_kinds:
        er_dia_result = run_coulomb_dia(
            wf,
            adia_densities,
            adia_ens,
            edmiston_ruedenberg_diabatization_jacobi,
            coulomb_tensor_fn,
        )
        # TODO: the six lines below could also be handled inside a function,
        # that is reused for all diabatization approaches.
        postprocess_dia_wrapped("er", er_dia_result)

    # Edmiston-Ruedenberg-ETA diabatization
    if DiaKind.EDMISTON_RUEDENBERG_ETA in dia_kinds:
        ereta_dia_result = run_coulomb_dia(
            wf,
            adia_densities,
            adia_ens,
            edmiston_ruedenberg_eta_diabatization,
            coulomb_tensor_fn,
            **dia_kwargs,
        )
        postprocess_dia_wrapped("ereta", ereta_dia_result, add_attrs=dia_kwargs)

    # Boys diabatization
    if DiaKind.BOYS in dia_kinds:
        boys_dia_result = run_boys_dia(dip_moms_2d, adia_ens)
        postprocess_dia_wrapped("boys", boys_dia_result)

    if CubeKind.DIA_DA in cube_kinds and 0 in states:
        cube_kinds &= ~CubeKind.DIA_DA
        logger.info(
            "Ground state involved! Disabled generation of diabatic "
            "detachment/attachment densities!"
        )

    # Loop over all diabatization that were carried out and calculate diabatic
    # properties as requested.
    for prefix, dia_result in all_dia_results.items():
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
                out_dir,
                reorder_coeffs=c2s_coeffs,
                **grid_kwargs,
            )
            cube_fns[f"{prefix.upper()}_DIA_SD"] = dia_sd_cube_fns

    return all_dia_results, cube_fns


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
    logger.info(f"Read wavefunction: {wf}")

    cis_fn = base.with_suffix(".cis")
    log_fn = base.with_suffix(".out")
    if not log_fn.exists():
        log_fn = base.with_suffix(".log")

    all_ens = parse_orca_all_energies(log_fn, do_tddft=True, triplets=triplets)
    logger.info("Parsed energies from ORCA log file")

    Xa, Ya, Xb, Yb = parse_orca_cis(
        cis_fn, restricted_same_ab=True, triplets_only=triplets
    )
    logger.info("Parsed ORCA .cis file")
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


def run(
    orca_outs: list[str],
    states: list[int],
    triplets: bool,
    dia_kinds: DiaKind,
    cube_kinds: CubeKind,
    ovlp: bool,
    grid_kwargs: dict,
    dia_kwargs: dict,
    out_dir: Path,
):
    nouts = len(orca_outs)
    nstates = len(states)

    all_dia_results = list()
    all_coords = list()
    dia_inp_prev = None
    all_states = np.zeros((nouts, nstates))
    for i, base in enumerate(orca_outs):
        base_name = base.stem
        # TODO: allow other codes
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
        if i > 0:
            logger.info(f"States at previous step were {dia_inp.states}.")
        # Update states by calculating overlaps if requested
        if ovlp:
            # TODO: report previous & new state before/after overlaps
            dia_inp.update_states(dia_inp_prev)
        logger.info(f"States for {base} are {dia_inp.states}")
        all_states[i] = states

        dia_results, cube_fns = run_dia(
            dia_inp,
            dia_kinds=dia_kinds,
            cube_kinds=cube_kinds,
            dia_kwargs=dia_kwargs,
            grid_kwargs=grid_kwargs,
            out_dir=out_dir,
        )
        all_dia_results.append(dia_results)
        all_coords.append(wf.coords)
        dia_inp_prev = dia_inp
        cb_fn = chimerabatch.write_inp(
            cube_fns, wf.atoms, wf.coords, base_name, out_dir
        )
        logger.info(f"Created script for chimerabatch.py in '{cb_fn}'.")

    np.save("all_states.npy", all_states)

    if nouts > 1:
        cds = get_coords_diffs(all_coords, align=True)
        for key, dia_result in dia_results.items():
            plot_dia_results(cds, dia_result, states, key)


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
    parser.add_argument(
        "-T",
        "--temperature",
        type=float,
        help="Temperature in K. Only required for EDMISTON_RUEDENBERG_ETA.",
    )
    parser.add_argument(
        "-C",
        "--pekar",
        type=float,
        help="Pekar factor. Only required for EDMISTON_RUEDENBERG_ETA.",
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
    parser.add_argument(
        "--out-dir",
        default=".",
    )

    # return parser.parse_args(args)
    args = parser.parse_args(args)

    if (DiaKind.EDMISTON_RUEDENBERG_ETA in args.dia) and (
        args.temperature is None or args.pekar is None
    ):
        parser.error(
            "When EDMISTON_RUEDENBERG_ETA is selected temperature (-T) "
            "and Pekar factor (-C) must also be given!"
        )
    return args


def run_cli():
    args = parse_args(sys.argv[1:])

    orca_outs = list(map(Path, args.orca_outs))
    states = args.states
    triplets = args.triplets
    dia_kinds = functools.reduce(operator.or_, args.dia)
    cube_kinds = functools.reduce(operator.or_, args.cube)
    grid_points = args.grid_points
    ovlp = args.ovlp
    out_dir = Path(args.out_dir)
    temperature = args.temperature
    pekar = args.pekar

    grid_kwargs = {
        "num": grid_points,
    }

    dia_kwargs = {
        "temperature": temperature,
        "pekar": pekar,
    }

    grid_kwargs = {
        "num": grid_points,
    }

    dia_kwargs = {
        "temperature": temperature,
        "pekar": pekar,
    }

    grid_points = args.grid_points
    out_dir = Path(args.out_dir)

    return run(
        orca_outs,
        states,
        triplets,
        dia_kinds,
        cube_kinds,
        ovlp,
        grid_kwargs,
        dia_kwargs,
        out_dir,
    )


if __name__ == "__main__":
    run_cli()
