# [1] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b01148
#     Plasser, 2016
# [2] https://doi.org/10.1002/jcc.25800
#     Garcia, Campetella, 2019
# [3] https://doi.org/10.1021/acs.jctc.0c00295
#     First Principles Nonadiabatic Excited-State Molecular Dynamics in NWChem
#     Song, Fischer, Zhang, Cramer, Mukamel, Govind, Tretiak, 2020

from collections import namedtuple
from pathlib import Path, PosixPath
import shutil
import tempfile

import h5py
import numpy as np
from scipy.optimize import linear_sum_assignment

from pysisyphus import logger
from pysisyphus.calculators.Calculator import Calculator

from pysisyphus.calculators.WFOWrapper import WFOWrapper
from pysisyphus.config import get_cmd
from pysisyphus.helpers_pure import describe
from pysisyphus.io.hdf5 import get_h5_group
from pysisyphus.wrapper.mwfn import make_cdd, get_mwfn_exc_str
from pysisyphus.wrapper.jmol import render_cdd_cube as render_cdd_cube_jmol
from pysisyphus.wavefunction.excited_states import (
    # nto_overlaps,
    norm_ci_coeffs,
    # nto_org_overlaps,
    tden_overlaps,
    top_differences,
)


NTOs = namedtuple("NTOs", "ntos lambdas")


def get_data_model(
    exc_state_num,
    occ_a,
    virt_a,
    occ_b,
    virt_b,
    ovlp_type,
    atoms,
    max_cycles,
):
    state_num = exc_state_num + 1  # including GS
    _1d = (max_cycles,)
    ovlp_state_num = state_num if ovlp_type == "wfow" else exc_state_num

    assert (occ_a + virt_a) == (occ_b + virt_b)
    nbf = occ_a + virt_a
    data_model = {
        "Ca": (max_cycles, nbf, nbf),
        "Cb": (max_cycles, nbf, nbf),
        "Xa": (max_cycles, exc_state_num, occ_a, virt_a),
        "Ya": (max_cycles, exc_state_num, occ_a, virt_a),
        "Xb": (max_cycles, exc_state_num, occ_b, virt_b),
        "Yb": (max_cycles, exc_state_num, occ_b, virt_b),
        "coords": (max_cycles, len(atoms) * 3),
        "all_energies": (
            max_cycles,
            state_num,
        ),
        "calculated_roots": _1d,
        "roots": _1d,
        "root_flips": _1d,
        "row_inds": _1d,
        "ref_cycles": _1d,
        "ref_roots": _1d,
        "overlap_matrices": (max_cycles, ovlp_state_num, ovlp_state_num),
        "cdd_cubes": _1d,
        "cdd_imgs": _1d,
    }
    return data_model


class GroundStateContext:
    def __init__(self, calc):
        self.calc = calc

    def __enter__(self):
        try:
            self.root_bak = self.calc.root
            self.calc.root = None
            self.track_bak = self.calc.track
            self.calc.track = False
        except AttributeError:
            pass
        try:
            self.wavefunction_dump_bak = self.calc.wavefunction_dump
            self.calc.wavefunction_dump = True
        except AttributeError:
            pass

    def __exit__(self, exc_type, exc_value, exc_trackback):
        try:
            self.calc.root = self.root_bak
            self.calc.track = self.track_bak
        except AttributeError:
            pass
        try:
            self.calc.wavefunction_dump = self.wavefunction_dump_bak
        except AttributeError:
            pass


def get_tden_overlaps(Ca1, Cb1, Xa1, Ya1, Xb1, Yb1, Ca2, Cb2, Xa2, Ya2, Xb2, Yb2, S_AO):
    XpYa1 = Xa1 + Ya1
    XpYa2 = Xa2 + Ya2
    ovlp_mat_a = tden_overlaps(Ca1, XpYa1, Ca2, XpYa2, S_AO)

    XpYb1 = Xb1 + Yb1
    XpYb2 = Xb2 + Yb2
    ovlp_mat_b = tden_overlaps(Cb1, XpYb1, Cb2, XpYb2, S_AO)

    ovlp_mat = ovlp_mat_a + ovlp_mat_b
    return ovlp_mat


def get_top_differences(
    Ca1, Cb1, Xa1, Ya1, Xb1, Yb1, Ca2, Cb2, Xa2, Ya2, Xb2, Yb2, S_AO
):
    """Transition orbital projection."""
    S_MO_a = Ca1.T @ S_AO @ Ca2
    S_MO_b = Cb1.T @ S_AO @ Cb2

    alpha_diffs = top_differences(Xa1, Ya1, Xa2, Ya2, S_MO_a)
    beta_diffs = top_differences(Xb1, Yb1, Xb2, Yb2, S_MO_b)
    diffs = alpha_diffs + beta_diffs
    return diffs


def get_ovlp_mat(
    Ca1: np.ndarray,
    Cb1: np.ndarray,
    Xa1: np.ndarray,
    Ya1: np.ndarray,
    Xb1: np.ndarray,
    Yb1: np.ndarray,
    Ca2: np.ndarray,
    Cb2: np.ndarray,
    Xa2: np.ndarray,
    Ya2: np.ndarray,
    Xb2: np.ndarray,
    Yb2: np.ndarray,
    S_AO: np.ndarray,
    ovlp_type: str,
):
    assert ovlp_type in ("tden", "top")

    if ovlp_type == "tden":
        ovlp_mat = get_tden_overlaps(
            Ca1, Cb1, Xa1, Ya1, Xb1, Yb1, Ca2, Cb2, Xa2, Ya2, Xb2, Yb2, S_AO
        )
    # elif ovlp_type == "wf":
    # raise Exception("wf-overlaps are not yet implemented!")
    # elif ovlp_type == "nto":
    # raise Exception("nto-overlaps are not yet implemented!")
    # elif ovlp_type == "nto_org":
    # raise Exception("nto_org-overlaps are not yet implemented!")
    elif ovlp_type == "top":
        top_rs = get_top_differences(
            Ca1, Cb1, Xa1, Ya1, Xb1, Yb1, Ca2, Cb2, Xa2, Ya2, Xb2, Yb2, S_AO
        )
        # Convert differences to some kind of pseudo-overlap matrix
        ovlp_mat = 1.0 - top_rs
    return ovlp_mat


def root2_from_ovlp_mat_and_root1(ovlp_mat, root1, min_cost=True):
    assert root1 > 0, "This function does not yet handle root1 == 0!"
    # TODO: update this function to work with wf-overlaps.
    # Two ideas: force the user to already provide a row_ind, e.g. calculate root1 - 1
    # outside this function, or add a flag that indicates that the ovlp_mat includes
    # the ground state which would be root 0. In such cases we would not subtract
    # 1 from the provided root.
    row_ind = root1 - 1

    # Work with the absolute overlap matrix, as phase switches between wavefunction
    # can lead to negative overlaps.
    abs_ovlp_mat = np.abs(ovlp_mat)

    if min_cost:
        # Match all excited state of the current and the reference step to make the
        # assignment more reasonable and avoid any double assignments.
        _, col_inds = linear_sum_assignment(-abs_ovlp_mat)
        ref_root_row = ovlp_mat[row_ind]
        root2 = col_inds[row_ind]
    else:
        # Match the best root row wise/just pick the root w/ the highest overlap.
        ref_root_row = abs_ovlp_mat[row_ind]
        root2 = ref_root_row.argmax()
    root2 += 1
    return root2


def track_root(
    Ca1,
    Cb1,
    Xa1,
    Ya1,
    Xb1,
    Yb1,
    Ca2,
    Cb2,
    Xa2,
    Ya2,
    Xb2,
    Yb2,
    S_AO,
    ovlp_type,
    root1,
    min_cost=True,
):
    ovlp_mat = get_ovlp_mat(
        Ca1, Cb1, Xa1, Ya1, Xb1, Yb1, Ca2, Cb2, Xa2, Ya2, Xb2, Yb2, S_AO, ovlp_type
    )
    return root2_from_ovlp_mat_and_root1(ovlp_mat, root1, min_cost=min_cost)


def args_from_calc(calc):
    Ca = calc.Ca_list[-1]
    Cb = calc.Cb_list[-1]
    Xa = calc.Xa_list[-1]
    Ya = calc.Ya_list[-1]
    Xb = calc.Xb_list[-1]
    Yb = calc.Yb_list[-1]
    return Ca, Cb, Xa, Ya, Xb, Yb


def track_root_between_ovlp_cals(calc1, calc2, **kwargs):
    ovlp_type = calc1.ovlp_type
    S_AO = calc1.get_sao_from_mo_coeffs(calc1.Ca_list[-1])
    return track_root(
        *args_from_calc(calc1),
        *args_from_calc(calc2),
        S_AO,
        ovlp_type,
        calc1.root,
        **kwargs,
    )


class OverlapCalculator(Calculator):
    OVLP_TYPE_VERBOSE = {
        "wf": "wavefunction overlap",
        "tden": "transition density matrix overlap",
        "nto": "natural transition orbital overlap",
        # As described in 10.1002/jcc.25800
        "nto_org": "original natural transition orbital overlap",
        "top": "transition orbital pair overlap",
    }
    VALID_KEYS = [
        k for k in OVLP_TYPE_VERBOSE.keys()
    ]  # lgtm [py/non-iterable-in-for-loop]
    VALID_CDDS = (None, "calc", "render")
    H5_MAP = {
        "Ca": "Ca_list",
        "Cb": "Cb_list",
        "Xa": "Xa_list",
        "Ya": "Ya_list",
        "Xb": "Xb_list",
        "Yb": "Yb_list",
        "coords": "coords_list",
        "all_energies": "all_energies_list",
        "roots": "roots_list",
        "ref_roots": "reference_roots",
    }

    def __init__(
        self,
        *args,
        root=None,
        nroots=None,
        track=False,
        ovlp_type="tden",
        double_mol=False,
        ovlp_with="previous",
        # TODO: reenable XY for wfoverlap?!
        # XY="X+Y",
        adapt_args=(0.5, 0.3, 0.6),
        # use_ntos=4,
        # pr_nto=False,
        # nto_thresh=0.3,
        cdds=None,
        orient="",
        dump_fn="overlap_data.h5",
        h5_dump=False,
        ncore=0,
        conf_thresh=1e-3,
        mos_ref="cur",
        mos_renorm: bool = True,
        min_cost: bool = False,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.root = root
        self.nroots = nroots
        self.track = track

        # TODO: enable this, when all calculators implement self.root & self.nroots
        # if self.track:
        # assert self.root <= self.nroots, "'root' must be smaller " "than 'nroots'!"
        self.ovlp_type = ovlp_type
        assert (
            self.ovlp_type in self.OVLP_TYPE_VERBOSE.keys()
        ), f"Valid overlap types are {self.VALID_KEYS}"
        self.double_mol = double_mol
        assert ovlp_with in ("previous", "first", "adapt")
        self.ovlp_with = ovlp_with
        assert (self.ovlp_type, self.ovlp_with) != (
            "top",
            "adapat",
        ), "ovlp_type: top and ovlp_with: adapat are not yet compatible"
        self.adapt_args = np.abs(adapt_args, dtype=float)
        self.adpt_thresh, self.adpt_min, self.adpt_max = self.adapt_args
        # self.use_ntos = use_ntos
        # self.pr_nto = pr_nto
        # self.nto_thresh = nto_thresh
        self.cdds = cdds
        # When calculation/rendering of charge density differences (CDDs) is
        # requested check fore the required programs (Multiwfn/Jmol). If they're
        # not available, we fallback to a more sensible command and print a warning.
        msg = (
            "'cdds: {0}' requested, but {1} was not found! "
            "Falling back to 'cdds: {2}'!\nConsider defining the {1} "
            "command in '.pysisyphusrc'."
        )

        jmol_cmd = get_cmd("jmol")
        mwfn_cmd = get_cmd("mwfn")
        if (self.cdds == "render") and not jmol_cmd:
            logger.warning(msg.format(self.cdds, "Jmol", "calc"))
            self.cdds = "calc"
        if (self.cdds in ("calc", "render")) and not mwfn_cmd:
            logger.warning(msg.format(self.cdds, "Multiwfn", None))
            self.cdds = None
        self.log(f"cdds: {self.cdds}, jmol={jmol_cmd}, mwfn={mwfn_cmd}")
        assert self.cdds in self.VALID_CDDS
        self.orient = orient
        self.dump_fn = self.out_dir / dump_fn
        self.h5_dump = h5_dump
        self.ncore = int(ncore)
        self.conf_thresh = float(conf_thresh)
        self.mos_ref = mos_ref
        assert self.mos_ref in ("cur", "ref")
        self.mos_renorm = bool(mos_renorm)
        self.min_cost = bool(min_cost)

        assert self.ncore >= 0, f"{ncore=} must be a positive number!"

        # MO-coefficients; MOs are expected to be in rows.
        self.Ca_list = list()
        self.Cb_list = list()
        # Transition density matrices.
        self.Xa_list = list()
        self.Ya_list = list()
        self.Xb_list = list()
        self.Yb_list = list()

        self.wfow = None
        self.nto_list = list()
        self.coords_list = list()
        # This list will hold the root indices at the beginning of the cycle
        # before any overlap calculation.
        self.calculated_roots = list()
        # This list will hold the (potentially) updated root after an overlap
        # calculation and it may differ from the value stored in
        # self.calculated_roots.
        self.roots_list = list()
        # Roots at the reference states that are used for comparison
        self.reference_roots = list()
        self.cdd_cubes = list()
        self.cdd_imgs = list()
        self.all_energies_list = list()
        # Why did is there already False in the list? Probably related
        # to plotting...
        self.root_flips = [
            False,
        ]
        self.first_root = None
        self.overlap_matrices = list()
        self.row_inds = list()
        # The first overlap calculation can be done in cycle 1, and then
        # we compare cycle 1 to cycle 0, regardless of the ovlp_with.
        self.ref_cycle = 0
        self.ref_cycles = list()
        self.atoms = None

        if track:
            self.log(
                "Tracking excited states with "
                f"{self.OVLP_TYPE_VERBOSE[ovlp_type]}s "
                f"between the current and the {self.ovlp_with} geometry."
            )
            if self.ovlp_with == "adapt":
                self.log(f"Adapt args: {self.adapt_args}")

        self.h5_fn = self.out_dir / "ovlp_data.h5"
        self.h5_group_name = self.name
        # We can't initialize the HDF5 group as we don't know the shape of
        # atoms/coords yet. So we wait until after the first calculation.
        self.h5_cycles = 50

        self._data_model = None
        self._h5_initialized = False

    def get_h5_group(self):
        if not self._h5_initialized:
            reset = True
            self._h5_initialized = True
        else:
            reset = False
        h5_group = get_h5_group(
            self.h5_fn, self.h5_group_name, self.data_model, reset=reset
        )
        return h5_group

    @property
    def data_model(self):
        if self._data_model is None:
            max_cycles = self.h5_cycles
            # exc_state_num, occ_mo_num, virt_mo_num = self.ci_coeff_list[0].shape
            exc_state_num, occ_a, virt_a = self.Xa_list[0].shape
            exc_state_num, occ_b, virt_b = self.Xb_list[0].shape
            self._data_model = get_data_model(
                exc_state_num,
                occ_a,
                virt_a,
                occ_b,
                virt_b,
                self.ovlp_type,
                self.atoms,
                max_cycles,
            )
        return self._data_model

    @property
    def stored_calculations(self):
        assert (
            len(self.Ca_list)
            == len(self.Cb_list)
            == len(self.Xa_list)
            == len(self.Ya_list)
            == len(self.Xb_list)
            == len(self.Yb_list)
        )
        return len(self.Xa_list)

    def clear_stored_calculations(self):
        for lst in (
            self.nto_list,
            self.coords_list,
            self.calculated_roots,
            self.roots_list,
            self.all_energies_list,
            self.Ca_list,
            self.Cb_list,
            self.Xa_list,
            self.Ya_list,
            self.Xb_list,
            self.Yb_list,
        ):
            lst.clear()

    def get_ci_coeffs_for(self, ind):
        return [
            lst[ind] for lst in (self.Xa_list, self.Ya_list, self.Xb_list, self.Yb_list)
        ]

    def prepare_overlap_data(self, path) -> tuple[
        np.ndarray,  # Ca
        np.ndarray,  # Xa
        np.ndarray,  # Ya
        np.ndarray,  # Cb
        np.ndarray,  # Xb
        np.ndarray,  # Yb
        np.ndarray,  # all_energies
    ]:
        """This method has to implement the calculator specific parsing of
        MO-coefficients, CI-coefficients and energies.
        Should return a filename pointing to TURBOMOLE
        like mos, a MO coefficient array and a CI coefficient array."""
        raise Exception("Implement me!")

    def store_overlap_data(self, atoms, coords, path=None, overlap_data=None):
        if self.atoms is None:
            self.atoms = atoms

        if overlap_data is None:
            overlap_data = self.prepare_overlap_data(path)

        if len(overlap_data) == 4:
            Ca, Xa, Ya, all_ens = overlap_data
            # Use same data for beta part
            Xb = Xa.copy()
            Yb = Ya.copy()
            Cb = Ca.copy()
        elif len(overlap_data) == 7:
            Ca, Xa, Ya, Cb, Xb, Yb, all_ens = overlap_data
        else:
            raise Exception("Expecting either 4 or 7 items in overlap_data!")

        assert Ca.ndim == Cb.ndim == 2
        assert all([mat.ndim == 3 for mat in (Xa, Ya, Xb, Yb)])

        Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)
        self.Ca_list.append(Ca.copy())
        self.Cb_list.append(Cb.copy())
        self.Xa_list.append(Xa)
        self.Ya_list.append(Ya)
        self.Xb_list.append(Xb)
        self.Yb_list.append(Yb)

        # If WF-overlaps were requested we initialize the WFOWrapper object
        if (self.ovlp_type == "wf") and (self.wfow is None):
            _, occa, virta = Xa.shape
            _, occb, virtb = Xb.shape
            self.set_wfow(occa, virta, occb, virtb)

        if self.first_root is None:
            self.first_root = self.root
            self.log(f"Set first root to {self.first_root}.")

        self.coords_list.append(coords)
        self.calculated_roots.append(self.root)
        # We can't calculate any overlaps in the first cycle, so we can't
        # compute a new root value. So we store the same value as for
        # calculated_roots.
        if self.stored_calculations < 2:
            self.roots_list.append(self.root)
        self.all_energies_list.append(all_ens)

        # Also store NTOs if requested
        # TODO: handle this differently
        # if self.ovlp_type in ("nto", "nto_org"):
        # self.set_ntos(mo_coeffs, ci_coeffs)

    def get_indices(self, indices=None):
        """
        A new root is determined by selecting the overlap matrix row
        corresponding to the reference root and checking for the root
        with the highest overlap (at the current geometry).

        The overlap matrix is usually formed by a double loop like:

        overlap_matrix = np.empty((ref_states, cur_states))
        for i, ref_state in enumerate(ref_states):
            for j, cur_state in enumerate(cur_states):
                overlap_matrix[i, j] = make_overlap(ref_state, cur_state)

        So the reference states run along the rows. Thats why the ref_state index
        comes first in the 'indices' tuple.
        """

        if indices is None:
            # By default we compare a reference cycle with the current (last)
            # cycle, so the second index is -1.
            ref, cur = self.ref_cycle, -1
        else:
            assert len(indices) == 2
            ref, cur = [int(i) for i in indices]
        return (ref, cur)

    @staticmethod
    def get_mo_norms(C, S_AO):
        return np.diag(C.T @ S_AO @ C)

    @staticmethod
    def renorm_mos(C, S_AO):
        norms = OverlapCalculator.get_mo_norms(C, S_AO)
        sqrts = np.sqrt(norms)
        return C / sqrts[None, :]

    def get_ref_mos(self, C_ref, C_cur):
        return {
            "ref": C_ref,
            "cur": C_cur,
        }[self.mos_ref]

    def get_orbital_matrices(self, indices=None, S_AO=None):
        """Return MO coefficents and AO overlaps for the given indices.

        If not provided, a AO overlap matrix is constructed from one of
        the MO coefficient matrices (controlled by self.mos_ref). Also,
        if requested one of the two MO coefficient matrices is re-normalized.
        """

        ref, cur = self.get_indices(indices)
        Ca_ref = self.Ca_list[ref].copy()
        Ca_cur = self.Ca_list[cur].copy()
        Cb_ref = self.Cb_list[ref].copy()
        Cb_cur = self.Cb_list[cur].copy()

        # Reconstruct from alpha MO coefficients
        if reconstruct_S_AO := (S_AO is None):
            C_S_AO = Ca_cur if (self.mos_ref == "cur") else Ca_ref
            self.log(f"Reconstructed S_AO from '{self.mos_ref}' MO coefficients.")
            S_AO = self.get_sao_from_mo_coeffs(C_S_AO)
            self.log(f"max(abs(S_AO))={np.abs(S_AO).max():.6f}")

        Cs = [[Ca_ref, Cb_ref], [Ca_cur, Cb_cur]]
        # Only renormalize if requested and we reconstructed the AO overlap matrix.
        if self.mos_renorm and reconstruct_S_AO:
            # If S_AO was reconstructed from "cur" MOs, then "ref" MOs won't be
            # normalized anymore and vice versa.
            renorm_ind = 0 if (self.mos_ref == "cur") else 1
            Cs[renorm_ind] = [self.renorm_mos(C, S_AO) for C in Cs[renorm_ind]]
            self.log(f"Renormalized '{('ref', 'cur')[renorm_ind]}' MO coefficients.")
        elif self.mos_renorm and (not reconstruct_S_AO):
            self.log("Skipped MO re-normalization as 'S_AO' was provided.")

        norms_0a = self.get_mo_norms(Cs[0][0], S_AO)
        norms_0b = self.get_mo_norms(Cs[0][1], S_AO)
        norms_1a = self.get_mo_norms(Cs[1][0], S_AO)
        norms_1b = self.get_mo_norms(Cs[1][1], S_AO)
        self.log(f"norm(MOs_0a): {describe(norms_0a)}")
        self.log(f"norm(MOs_0b): {describe(norms_0b)}")
        self.log(f"norm(MOs_1a): {describe(norms_1a)}")
        self.log(f"norm(MOs_1b): {describe(norms_1b)}")
        return *Cs[0], *Cs[1], S_AO

    @staticmethod
    def get_sao_from_mo_coeffs(C):
        """Recover AO overlaps from given MO coefficients.

        For MOs in the columns of mo_coeffs:

            S_AO = C⁻¹^T C⁻¹
            S_AO C = C⁻¹^T
            (S_AO C)^T = C⁻¹
            C^T S_AO^T = C⁻¹
            C^T S_AO C = I
        """
        C_inv = np.linalg.pinv(C, rcond=1e-8)
        S_AO = C_inv.T @ C_inv
        return S_AO

    def get_wf_overlaps(self, indices=None, S_AO=None):
        Ca_ref, Cb_ref, Ca_cur, Cb_cur, S_AO = self.get_orbital_matrices(indices, S_AO)
        C_ref = np.concatenate((Ca_ref, Cb_ref), axis=1)
        C_cur = np.concatenate((Ca_cur, Cb_cur), axis=1)

        ref, cur = self.get_indices(indices)
        # Reference step
        Xa_ref, Ya_ref, Xb_ref, Yb_ref = self.get_ci_coeffs_for(ref)
        # Current step
        Xa_cur, Ya_cur, Xb_cur, Yb_cur = self.get_ci_coeffs_for(cur)
        ref_cycle = (C_ref, Xa_ref + Ya_ref, Xb_ref + Yb_ref)
        cur_cycle = (C_cur, Xa_cur + Ya_cur, Xb_cur + Yb_cur)
        return self.wfow.wf_overlap(*ref_cycle, *cur_cycle, ao_ovlp=S_AO)

    """
    def wf_overlaps(self, mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2, S_AO=None):
        cycle1 = (mo_coeffs1, ci_coeffs1)
        cycle2 = (mo_coeffs2, ci_coeffs2)
        overlaps = self.wfow.wf_overlap(cycle1, cycle2, S_AO=S_AO)
        return overlaps

    def wf_overlap_with_calculator(self, calc, S_AO=None):
        mo_coeffs1 = self.mo_coeff_list[-1]
        ci_coeffs1 = self.ci_coeff_list[-1]
        mo_coeffs2 = calc.mo_coeff_list[-1]
        ci_coeffs2 = calc.ci_coeff_list[-1]
        overlaps = self.wf_overlaps(
            mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2, S_AO=S_AO
        )
        return overlaps
    """

    def get_tden_overlaps(self, indices=None, S_AO=None):
        Ca_ref, Cb_ref, Ca_cur, Cb_cur, S_AO = self.get_orbital_matrices(indices, S_AO)

        ref, cur = self.get_indices(indices)
        # Reference step
        Xa_ref, Ya_ref, Xb_ref, Yb_ref = self.get_ci_coeffs_for(ref)
        # Current step
        Xa_cur, Ya_cur, Xb_cur, Yb_cur = self.get_ci_coeffs_for(cur)

        return get_tden_overlaps(
            Ca_ref,
            Cb_ref,
            Xa_ref,
            Ya_ref,
            Xb_ref,
            Yb_ref,
            Ca_cur,
            Cb_cur,
            Xa_cur,
            Ya_cur,
            Xb_cur,
            Yb_cur,
            S_AO,
        )

    """
    def calculate_state_ntos(self, state_ci_coeffs, mos):
        # TODO: don't renorm; respect self.XY choice 
        normed = state_ci_coeffs / np.linalg.norm(state_ci_coeffs)
        # u, s, vh = np.linalg.svd(state_ci_coeffs)
        u, s, vh = np.linalg.svd(normed)
        lambdas = s ** 2
        self.log("Normalized transition density vector to 1.")
        self.log(f"Sum(lambdas)={np.sum(lambdas):.4f}")
        lambdas_str = np.array2string(lambdas[:3], precision=4, suppress_small=True)
        self.log(f"First three lambdas: {lambdas_str}")

        occ_mo_num = state_ci_coeffs.shape[0]
        occ_mos = mos[:occ_mo_num]
        vir_mos = mos[occ_mo_num:]
        occ_ntos = occ_mos.T.dot(u)
        vir_ntos = vir_mos.T.dot(vh)
        return occ_ntos, vir_ntos, lambdas

    def get_nto_overlaps(self, indices=None, S_AO=None, org=False):
        ref, cur = self.get_indices(indices)

        if S_AO is None:
            S_AO = self.get_sao_from_mo_coeffs_and_dump(
                self.get_ref_mos(self.mo_coeff_list[ref], self.mo_coeff_list[cur])
            )

        ntos_1 = self.nto_list[ref]
        ntos_2 = self.nto_list[cur]
        if org:
            overlaps = nto_org_overlaps(
                ntos_1, ntos_2, S_AO, nto_thresh=self.nto_thresh
            )
        else:
            overlaps = nto_overlaps(ntos_1, ntos_2, S_AO)
        return overlaps
    """

    """
    def set_ntos(self, mo_coeffs, ci_coeffs):
        roots = ci_coeffs.shape[0]
        ntos_for_cycle = list()
        for root in range(roots):
            sn_ci_coeffs = ci_coeffs[root]
            self.log(f"Calculating NTOs for root {root+1}")
            occ_ntos, vir_ntos, lambdas = self.calculate_state_ntos(
                sn_ci_coeffs,
                mo_coeffs,
            )
            pr_nto = lambdas.sum() ** 2 / (lambdas ** 2).sum()
            if self.pr_nto:
                use_ntos = int(np.round(pr_nto))
                self.log(f"PR_NTO={pr_nto:.2f}")
            else:
                use_ntos = self.use_ntos
            self.log(f"Using {use_ntos} NTOS")
            ovlp_occ_ntos = occ_ntos.T[:use_ntos]
            ovlp_vir_ntos = vir_ntos.T[:use_ntos]
            ovlp_lambdas = lambdas[:use_ntos]
            ovlp_lambdas = np.concatenate((ovlp_lambdas, ovlp_lambdas))
            ovlp_ntos = np.concatenate((ovlp_occ_ntos, ovlp_vir_ntos), axis=0)
            ntos = NTOs(ntos=ovlp_ntos, lambdas=ovlp_lambdas)
            ntos_for_cycle.append(ntos)
        self.nto_list.append(ntos_for_cycle)
    """

    def get_top_differences(self, indices=None, S_AO=None):
        """Transition orbital projection."""
        Ca_ref, Cb_ref, Ca_cur, Cb_cur, S_AO = self.get_orbital_matrices(indices, S_AO)
        S_MO_a = Ca_ref.T @ S_AO @ Ca_cur
        S_MO_b = Cb_ref.T @ S_AO @ Cb_cur

        ref, cur = self.get_indices(indices)
        # Reference step
        Xa_ref, Ya_ref, Xb_ref, Yb_ref = self.get_ci_coeffs_for(ref)
        # Current step
        Xa_cur, Ya_cur, Xb_cur, Yb_cur = self.get_ci_coeffs_for(cur)

        alpha_diffs = top_differences(Xa_ref, Ya_ref, Xa_cur, Ya_cur, S_MO_a)
        beta_diffs = top_differences(Xb_ref, Yb_ref, Xb_cur, Yb_cur, S_MO_b)
        diffs = alpha_diffs + beta_diffs
        return diffs

    def dump_overlap_data(self):
        if self.h5_dump:
            h5_group = self.get_h5_group()

            h5_group.attrs["ovlp_type"] = self.ovlp_type
            h5_group.attrs["ovlp_with"] = self.ovlp_with
            h5_group.attrs["orient"] = self.orient
            h5_group.attrs["atoms"] = np.bytes_(self.atoms)

            for key, shape in self.data_model.items():
                try:
                    mapped_key = self.H5_MAP[key]
                except KeyError:
                    mapped_key = key
                value = getattr(self, mapped_key)
                # Skip this value if the underlying list is empty
                if not value:
                    continue
                cur_cycle = self.calc_counter
                cur_value = value[-1]
                # the CDD strings are currently not yet handled properly
                if type(cur_value) == PosixPath:
                    cur_value = str(cur_value)
                    continue
                if len(shape) > 1:
                    h5_group[key][cur_cycle, : len(cur_value)] = cur_value
                else:
                    h5_group[key][cur_cycle] = cur_value

        data_dict = {
            "Ca": np.array(self.Ca_list, dtype=float),
            "Cb": np.array(self.Cb_list, dtype=float),
            "Xa": np.array(self.Xa_list, dtype=float),
            "Ya": np.array(self.Ya_list, dtype=float),
            "Xb": np.array(self.Xb_list, dtype=float),
            "Yb": np.array(self.Yb_list, dtype=float),
            "coords": np.array(self.coords_list, dtype=float),
            "all_energies": np.array(self.all_energies_list, dtype=float),
        }
        if self.root:
            root_dict = {
                "calculated_roots": np.array(self.calculated_roots, dtype=int),
                "roots": np.array(self.roots_list, dtype=int),
                "root_flips": np.array(self.root_flips, dtype=bool),
                "overlap_matrices": np.array(self.overlap_matrices, dtype=float),
                "row_inds": np.array(self.row_inds, dtype=int),
                "ref_cycles": np.array(self.ref_cycles, dtype=int),
                "ref_roots": np.array(self.reference_roots, dtype=int),
            }
            data_dict.update(root_dict)

        if self.cdd_cubes:
            data_dict["cdd_cubes"] = np.array(self.cdd_cubes, dtype="S")
            if self.cdd_imgs:
                data_dict["cdd_imgs"] = np.array(self.cdd_imgs, dtype="S")

        with h5py.File(self.dump_fn, "w") as handle:
            for key, val in data_dict.items():
                if key in ("Ca", "Cb", "Xa", "Ya", "Xb", "Yb"):
                    add_kwargs = {
                        "compression": "gzip",
                        "compression_opts": 9,
                    }
                else:
                    add_kwargs = {}
                handle.create_dataset(name=key, dtype=val.dtype, data=val, **add_kwargs)
            handle.attrs["ovlp_type"] = self.ovlp_type
            handle.attrs["ovlp_with"] = self.ovlp_with
            handle.attrs["orient"] = self.orient
            handle.attrs["atoms"] = np.array(self.atoms, "S1")

    @staticmethod
    def from_overlap_data(h5_fn, set_wfow=False):
        calc = OverlapCalculator(track=True)

        root_info = False
        with h5py.File(h5_fn) as handle:
            try:
                ovlp_with = handle["ovlp_with"][()].decode()
                ovlp_type = handle["ovlp_type"][()].decode()
            except KeyError:
                ovlp_with = handle.attrs["ovlp_with"]
                ovlp_type = handle.attrs["ovlp_type"]
            all_energies = handle["all_energies"][:]
            try:
                ref_roots = handle["ref_roots"][:]
                roots = handle["roots"][:]
                calculated_roots = handle["calculated_roots"][:]
                root_info = True
            except KeyError:
                print(f"Couldn't find root information in '{h5_fn}'.")
            # TODO: set using H5 map?
            calc.Ca_list = handle["Ca"][:]
            calc.Cb_list = handle["Cb"][:]
            calc.Xa_list = handle["Xa"][:]
            calc.Ya_list = handle["Ya"][:]
            calc.Xb_list = handle["Xb"][:]
            calc.Yb_list = handle["Yb"][:]

        calc.ovlp_type = ovlp_type
        calc.ovlp_with = ovlp_with
        calc.all_energies_list = list(all_energies)
        if root_info:
            calc.roots_list = list(roots)
            calc.calculated_roots = list(calculated_roots)
            try:
                calc.first_root = ref_roots[0]
                calc.root = calc.first_root
            except IndexError:
                calc.root = roots[0]

        # TODO: set wfow
        # if (ovlp_type == "wf") or set_wfow:
        # calc.set_wfow(ci_coeffs[0])

        return calc

    def set_wfow(self, occa: int, virta: int, occb: int, virtb: int):
        assert occa + virta == occb + virtb
        try:
            wfow_mem = self.pal * self.mem
        except AttributeError:
            wfow_mem = 8000
        self.wfow = WFOWrapper(
            occa,
            virta,
            occb,
            virtb,
            calc_number=self.calc_number,
            wfow_mem=wfow_mem,
            ncore=self.ncore,
            conf_thresh=self.conf_thresh,
            out_dir=self.out_dir,
        )

    def track_root(self, ovlp_type=None):
        """Check if a root flip occured occured compared to the previous cycle
        by calculating the overlap matrix wrt. a reference cycle."""

        if ovlp_type is None:
            ovlp_type = self.ovlp_type

        # Nothing to compare to if only one calculation was done yet OR self.root
        # is set to None, e.g., there is nothing to compare to.
        # Nonetheless, dump the first cycle to HDF5.
        if (self.root is None) or (self.stored_calculations < 2):
            self.dump_overlap_data()
            self.log(
                "Skipping overlap calculation in the first cycle "
                "as there is nothing to compare to."
            )
            return False

        S_AO = None
        # We can only run a double molecule calculation if it is
        # implemented for the specific calculator.
        if self.double_mol and hasattr(self, "run_double_mol_calculation"):
            old, new = self.get_indices()
            two_coords = self.coords_list[old], self.coords_list[new]
            S_AO = self.run_double_mol_calculation(self.atoms, *two_coords)
        elif (self.double_mol is False) and (self.ovlp_type == "wf"):
            # TODO: respect ref_mos?!
            S_AO = self.get_sao_from_mo_coeffs(self.Ca_list[-1])
            self.log("Created S_AO to avoid its creation in WFOverlap.")

        self.log(f"Calculating '{self.ovlp_type}' overlaps.")
        if ovlp_type == "wf":
            overlap_mats = self.get_wf_overlaps(S_AO=S_AO)
            # Use SVD-orthogonalized overlap matrix
            overlaps = np.abs(overlap_mats[2])
        elif ovlp_type == "tden":
            overlaps = self.get_tden_overlaps(S_AO=S_AO)
        elif ovlp_type == "nto":
            raise Exception("nto-overlaps are not yet implemented!")
            overlaps = self.get_nto_overlaps(S_AO=S_AO)
        elif ovlp_type == "nto_org":
            raise Exception("nto_org-overlaps are not yet implemented!")
            overlaps = self.get_nto_overlaps(S_AO=S_AO, org=True)
        elif ovlp_type == "top":
            top_rs = self.get_top_differences(S_AO=S_AO)
            overlaps = 1.0 - top_rs
        else:
            raise Exception(
                "Invalid overlap type key! Use one of " + ", ".join(self.VALID_KEYS)
            )
        self.overlap_matrices.append(overlaps)
        overlaps = np.abs(overlaps)

        # In the end we are looking for a root flip compared to the
        # previous cycle.
        # This is done by comparing the excited states at the current cycle
        # to some older cycle (the first, the previous, or some cycle
        # in between), by determining the highest overlap in a given row
        # of the overlap matrix.
        ref_root = self.roots_list[self.ref_cycle]
        self.reference_roots.append(ref_root)
        # Row index in the overlap matrix. Depends on the root of the reference
        # cycle and corresponds to the old root.
        row_ind = ref_root - 1
        # With WFOverlaps the ground state is also present and the overlap
        # matrix has shape (N+1, N+1) instead of (N, N), with N being the
        # number of excited states.
        if self.ovlp_type == "wf":
            row_ind += 1
        self.row_inds.append(row_ind)
        self.ref_cycles.append(self.ref_cycle)
        self.log(f"Reference is cycle {self.ref_cycle}, root {ref_root}.")

        # As described in [3].
        # Code contributed by PT.
        if self.min_cost:
            # Match all excited state of the current and the reference step to make the
            # assignment more reasonable and avoid any double assignments.
            self.log("Assigned roots using Kuhn-Munkres algorithm.")
            _, col_inds = linear_sum_assignment(-overlaps)
            ref_root_row = overlaps[row_ind]
            new_root = col_inds[row_ind]
            if ref_root_row.argmax() != new_root:
                self.log(
                    "The newly assigned root is not root of highest overlap because "
                    "another row mapping to the same state had a higher overlap!"
                )
        else:
            # Match the best root row wise/just pick the root w/ the highest overlap.
            self.log(f"Analyzing row {row_ind} of the overlap matrix.")
            ref_root_row = overlaps[row_ind]
            new_root = ref_root_row.argmax()

        max_overlap = ref_root_row[new_root]
        if self.ovlp_type == "wf":
            new_root -= 1
        prev_root = self.root
        self.log(f"Root at previous cycle is {prev_root}.")
        self.root = new_root + 1
        ref_root_row_str = ", ".join(
            [f"{i}: {ov:.2%}" for i, ov in enumerate(ref_root_row)]
        )
        self.log(f"Overlaps: {ref_root_row_str}")
        root_flip = self.root != prev_root
        self.log(f"Highest overlap is {max_overlap:.2%}.")
        if not root_flip:
            self.log(f"Keeping current root {self.root}.")
        else:
            self.log(
                f"Root flip! New root is {self.root}. Root at previous "
                f"step was {prev_root}."
            )
        # Look for a new reference state if requested. We want to avoid
        # overlap matrices just after a root flip.
        if self.ovlp_with == "previous":
            self.ref_cycle += 1
        elif (self.ovlp_with == "adapt") and not root_flip:
            self.log("Checking wether the reference cycle has to be adapted.")
            sorted_inds = ref_root_row.argsort()
            sec_highest, highest = ref_root_row[sorted_inds[-2:]]
            ratio = sec_highest / highest
            self.log(
                f"Two highest overlaps: {sec_highest:.2%}, {highest:.2%}, "
                f"ratio={ratio:.4f}"
            )
            above_thresh = highest >= self.adpt_thresh
            self.log(
                f"Highest overlap is above threshold? (>= {self.adpt_thresh:.4f}): "
                f"{above_thresh}"
            )
            valid_ratio = self.adpt_min < ratio < self.adpt_max
            self.log(
                f"Ratio is valid? (between {self.adpt_min:.4f} and "
                f"{self.adpt_max:.4f}): {valid_ratio}"
            )
            """Only adapt the reference cycle when the overlaps are well
            behaved and the following two conditions are True:

            1.) The highest overlap is above the threshold.

            2.) The ratio value of (second highest)/(highest) is valid. A small
            value indicates cleary separated states and we probably
            don't have to update the reference cycle as the overlaps are still
            big enough.
            As the overlaps between two states become more similar the ratio
            approaches 1. This may occur in regions of state crossings and then
            we dont't want to update the reference cycle.
            """
            if above_thresh and valid_ratio:
                self.ref_cycle = len(self.calculated_roots) - 1

        if self.ref_cycle != self.ref_cycles[-1]:
            self.log(f"New reference cycle is {self.ref_cycle}.")
        else:
            self.log(f"Keeping old reference cycle {self.ref_cycle}.")

        self.root_flips.append(root_flip)
        self.roots_list.append(self.root)
        assert len(self.roots_list) == len(self.calculated_roots)

        if self.cdds:
            try:
                self.calc_cdd_cube(self.root)
            except Exception as err:
                print("CDD calculation by Multiwfn crashed. Disabling it!")
                self.log(err)
                self.cdds = None
            if self.cdds == "render":
                self.render_cdd_cube()

        self.dump_overlap_data()
        self.log("\n")

        # True if a root flip occured
        return root_flip

    def calc_cdd_cube(self, root, cycle=-1):
        # Check if Calculator provides an input file (.fchk/.molden) for Mwfn
        if not hasattr(self, "mwfn_wf"):
            self.log(
                "Calculator does not provide an input file for Multiwfn, "
                "as 'self.mwfn_wf' is not set! Skipping CDD cube generation!"
            )
        if cycle != -1:
            self.log("'cycle' argument to make_cdd_cube is currently ignored!")
        energies = self.all_energies_list[cycle]
        # As of Multiwfn 3.8 it seems that MWFN can't handle plain text input
        # with spin labels, so we only provide the alpha part..
        Xa, Ya, *_ = self.get_ci_coeffs_for(cycle)
        exc_str = get_mwfn_exc_str(energies, Xa=Xa, Ya=Ya)
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            exc_path = tmp_path / "exc_input"
            with open(exc_path, "w") as handle:
                handle.write(exc_str)
            cubes = make_cdd(self.mwfn_wf, root, str(exc_path), tmp_path)
            assert len(cubes) == 1
            cube = cubes[0]
            cube_fn = cubes[0].name
            new_cube_fn = self.make_fn(cube_fn)
            shutil.copy(cube, new_cube_fn)
            self.cdd_cubes.append(new_cube_fn)

    def render_cdd_cube(self):
        cdd_cube = self.cdd_cubes[-1]
        try:
            cdd_img = render_cdd_cube_jmol(cdd_cube, orient=self.orient)
            self.cdd_imgs.append(cdd_img)
        except:
            self.log("Something went wrong while rendering the CDD cube.")
