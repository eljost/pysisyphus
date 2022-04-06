# [1] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b01148
#     Plasser, 2016
# [2] https://doi.org/10.1002/jcc.25800
#     Garcia, Campetella, 2019

from collections import namedtuple
from pathlib import Path, PosixPath
import shutil
import tempfile

import h5py
import numpy as np

from pysisyphus import logger
from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.WFOWrapper import WFOWrapper
from pysisyphus.config import get_cmd
from pysisyphus.helpers_pure import describe
from pysisyphus.io.hdf5 import get_h5_group
from pysisyphus.wrapper.mwfn import make_cdd, get_mwfn_exc_str
from pysisyphus.wrapper.jmol import render_cdd_cube as render_cdd_cube_jmol


NTOs = namedtuple("NTOs", "ntos lambdas")


def get_data_model(
    exc_state_num, occ_mo_num, virt_mo_num, ovlp_type, atoms, max_cycles
):
    mo_num = occ_mo_num + virt_mo_num
    state_num = exc_state_num + 1  # including GS
    _1d = (max_cycles,)
    ovlp_state_num = state_num if ovlp_type == "wfow" else exc_state_num

    data_model = {
        "mo_coeffs": (max_cycles, mo_num, mo_num),
        "ci_coeffs": (max_cycles, exc_state_num, occ_mo_num, virt_mo_num),
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


class OverlapCalculator(Calculator):
    OVLP_TYPE_VERBOSE = {
        "wf": "wavefunction overlap",
        "tden": "transition density matrix overlap",
        "nto": "natural transition orbital overlap",
        # As described in 10.1002/jcc.25800
        "nto_org": "original natural transition orbital overlap",
    }
    VALID_KEYS = [
        k for k in OVLP_TYPE_VERBOSE.keys()
    ]  # lgtm [py/non-iterable-in-for-loop]
    VALID_CDDS = (None, "calc", "render")
    VALID_XY = ("X", "X+Y", "X-Y")
    H5_MAP = {
        "mo_coeffs": "mo_coeff_list",
        "ci_coeffs": "ci_coeff_list",
        "coords": "coords_list",
        "all_energies": "all_energies_list",
        "roots": "roots_list",
        "ref_roots": "reference_roots",
    }

    def __init__(
        self,
        *args,
        track=False,
        ovlp_type="tden",
        double_mol=False,
        ovlp_with="previous",
        XY="X+Y",
        adapt_args=(0.5, 0.3, 0.6),
        use_ntos=4,
        pr_nto=False,
        nto_thresh=0.3,
        cdds=None,
        orient="",
        dump_fn="overlap_data.h5",
        h5_dump=False,
        ncore=0,
        conf_thresh=1e-3,
        dyn_roots=0,
        mos_ref="cur",
        mos_renorm=True,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.track = track
        self.ovlp_type = ovlp_type
        assert (
            self.ovlp_type in self.OVLP_TYPE_VERBOSE.keys()
        ), f"Valid overlap types are {self.VALID_KEYS}"
        self.double_mol = double_mol
        assert ovlp_with in ("previous", "first", "adapt")
        self.ovlp_with = ovlp_with
        self.XY = XY
        assert self.XY in self.VALID_XY
        self.adapt_args = np.abs(adapt_args, dtype=float)
        self.adpt_thresh, self.adpt_min, self.adpt_max = self.adapt_args
        self.use_ntos = use_ntos
        self.pr_nto = pr_nto
        self.nto_thresh = nto_thresh
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
            logger.debug(msg.format(self.cdds, "Jmol", "calc"))
            self.cdds = "calc"
        if (self.cdds in ("calc", "render")) and not mwfn_cmd:
            logger.debug(msg.format(self.cdds, "Multiwfn", None))
            self.cdds = None
        self.log(f"cdds: {self.cdds}, jmol={jmol_cmd}, mwfn={mwfn_cmd}")
        assert self.cdds in self.VALID_CDDS
        self.orient = orient
        self.dump_fn = self.out_dir / dump_fn
        self.h5_dump = h5_dump
        self.ncore = int(ncore)
        self.conf_thresh = float(conf_thresh)
        self.dyn_roots = int(dyn_roots)
        if self.dyn_roots != 0:
            self.dyn_roots = 0
            self.log("dyn_roots = 0 is hardcoded right now")
        self.mos_ref = mos_ref
        assert self.mos_ref in ("cur", "ref")
        self.mos_renorm = bool(mos_renorm)

        assert self.ncore >= 0, "ncore must be a >= 0!"

        self.wfow = None
        self.mo_coeff_list = list()
        self.ci_coeff_list = list()
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
        self.root = None

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
            exc_state_num, occ_mo_num, virt_mo_num = self.ci_coeff_list[0].shape
            self._data_model = get_data_model(
                exc_state_num,
                occ_mo_num,
                virt_mo_num,
                self.ovlp_type,
                self.atoms,
                max_cycles,
            )
        return self._data_model

    @property
    def roots_number(self):
        return self.root + self.dyn_roots

    def blowup_ci_coeffs(self, ci_coeffs):
        states, occ, virt = ci_coeffs.shape
        full = np.zeros((states, occ, occ + virt))
        full[:, :, occ:] = ci_coeffs
        return full

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
    def get_mo_norms(mo_coeffs, ao_ovlp):
        """MOs are in rows."""
        # einsum-call is extremely slow
        # mo_norms = np.einsum("ki,kj,ij->k", mo_coeffs, mo_coeffs, ao_ovlp)
        mo_norms = np.diag(mo_coeffs.dot(ao_ovlp).dot(mo_coeffs.T))
        return mo_norms

    @staticmethod
    def renorm_mos(mo_coeffs, ao_ovlp):
        norms = OverlapCalculator.get_mo_norms(mo_coeffs, ao_ovlp)
        sqrts = np.sqrt(norms)
        return mo_coeffs / sqrts[:, None]

    def get_ref_mos(self, ref_mo_coeffs, cur_mo_coeffs):
        return {
            "ref": ref_mo_coeffs,
            "cur": cur_mo_coeffs,
        }[self.mos_ref]

    def get_orbital_matrices(self, indices=None, ao_ovlp=None):
        """Return MO coefficents and AO overlaps for the given indices.

        If not provided, a AO overlap matrix is constructed from one of
        the MO coefficient matrices (controlled by self.mos_ref). Also,
        if requested one of the two MO coefficient matrices is re-normalized.
        """

        ref, cur = self.get_indices(indices)
        ref_mo_coeffs = self.mo_coeff_list[ref].copy()
        cur_mo_coeffs = self.mo_coeff_list[cur].copy()

        ao_ovlp_reconstructed = ao_ovlp is None
        if ao_ovlp_reconstructed:
            sao_mo_coeffs = cur_mo_coeffs if (self.mos_ref == "cur") else ref_mo_coeffs
            self.log(f"Reconstructed S_AO from '{self.mos_ref}' MO coefficients.")
            ao_ovlp = self.get_sao_from_mo_coeffs(sao_mo_coeffs)
            self.log(f"max(abs(S_AO))={np.abs(ao_ovlp).max():.6f}")

        return_mos = [ref_mo_coeffs, cur_mo_coeffs]
        # Only renormalize if requested and we reconstructed the AO overlap matrix.
        if self.mos_renorm and ao_ovlp_reconstructed:
            # If S_AO was reconstructed from "cur" MOs, then "ref" MOs won't be
            # normalized anymore and vice versa.
            renorm_ind = 0 if (self.mos_ref == "cur") else 1
            to_renorm = return_mos[renorm_ind]
            # norms = self.get_mo_norms(to_renorm, ao_ovlp)
            return_mos[renorm_ind] = self.renorm_mos(to_renorm, ao_ovlp)
            self.log(f"Renormalized '{('ref', 'cur')[renorm_ind]}' MO coefficients.")
        elif self.mos_renorm and (not ao_ovlp_reconstructed):
            self.log("Skipped MO re-normalization as 'ao_ovlp' was provided.")

        # return *return_mos, ao_ovlp
        # The statement above is only valid in python>=3.8
        norms0 = self.get_mo_norms(return_mos[0], ao_ovlp)
        norms1 = self.get_mo_norms(return_mos[1], ao_ovlp)
        self.log(f"norm(MOs_0): {describe(norms0)}")
        self.log(f"norm(MOs_1): {describe(norms1)}")
        return return_mos[0], return_mos[1], ao_ovlp

    @staticmethod
    def get_sao_from_mo_coeffs(mo_coeffs):
        """Recover AO overlaps from given MO coefficients.

        For MOs in the columns of mo_coeffs:

            S_AO = C⁻¹^T C⁻¹
            S_AO C = C⁻¹^T
            (S_AO C)^T = C⁻¹
            C^T S_AO^T = C⁻¹
            C^T S_AO C = I

        Here, MOs are expected to be in rows of mo_coeffs, yielding

            C S_AO C^T = I
        """
        mo_coeffs_inv = np.linalg.pinv(mo_coeffs, rcond=1e-8)
        ao_ovlp = mo_coeffs_inv.dot(mo_coeffs_inv.T)
        return ao_ovlp

    def get_sao_from_mo_coeffs_and_dump(self, mo_coeffs):
        ao_ovlp = self.get_sao_from_mo_coeffs(mo_coeffs)
        ao_ovlp_fn = self.make_fn("ao_ovlp_rec")
        np.savetxt(ao_ovlp_fn, ao_ovlp)
        return ao_ovlp

    def get_wf_overlaps(self, indices=None, ao_ovlp=None):
        old, new = self.get_indices(indices)
        old_cycle = (self.mo_coeff_list[old], self.ci_coeff_list[old])
        new_cycle = (self.mo_coeff_list[new], self.ci_coeff_list[new])
        return self.wfow.wf_overlap(old_cycle, new_cycle, ao_ovlp)

    def wf_overlaps(self, mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2, ao_ovlp=None):
        cycle1 = (mo_coeffs1, ci_coeffs1)
        cycle2 = (mo_coeffs2, ci_coeffs2)
        overlaps = self.wfow.wf_overlap(cycle1, cycle2, ao_ovlp=ao_ovlp)
        return overlaps

    def wf_overlap_with_calculator(self, calc, ao_ovlp=None):
        mo_coeffs1 = self.mo_coeff_list[-1]
        ci_coeffs1 = self.ci_coeff_list[-1]
        mo_coeffs2 = calc.mo_coeff_list[-1]
        ci_coeffs2 = calc.ci_coeff_list[-1]
        overlaps = self.wf_overlaps(
            mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2, ao_ovlp=ao_ovlp
        )
        return overlaps

    def tden_overlaps(self, mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2, ao_ovlp):
        """
        Parameters
        ----------
        mo_coeffs1 : ndarray, shape (MOs, AOs)
            MO coefficient matrix. One row per MO, one column per basis
            function. Usually square.
        mo_coeffs2 : ndarray
            See mo_coeffs1.
        ci_coeffs1 : ndarray, shape(occ. MOs, MOs)
            CI-coefficient matrix.
        ci_coeffs2 : ndarray, shape(occ. MOs, MOs)
            See ci_coeffs1.
        ao_ovlp : ndarray, shape(AOs1, AOs2)
            Double molcule AO overlaps.
        """
        states1, occ, _ = ci_coeffs1.shape
        ci_full1 = self.blowup_ci_coeffs(ci_coeffs1)
        ci_full2 = self.blowup_ci_coeffs(ci_coeffs2)

        # MO overlaps
        S_MO = mo_coeffs1.dot(ao_ovlp).dot(mo_coeffs2.T)
        S_MO_occ = S_MO[:occ, :occ]

        overlaps = list()
        for state1 in ci_full1:
            precontr = S_MO_occ.dot(state1).dot(S_MO)
            for state2 in ci_full2:
                overlaps.append(np.sum(precontr * state2))
        overlaps = np.array(overlaps).reshape(states1, -1)

        return overlaps

    def get_tden_overlaps(self, indices=None, ao_ovlp=None):
        mo_coeffs_ref, mo_coeffs_cur, ao_ovlp = self.get_orbital_matrices(
            indices, ao_ovlp
        )

        ref, cur = self.get_indices(indices)
        ci_coeffs_ref = self.ci_coeff_list[ref]
        ci_coeffs_cur = self.ci_coeff_list[cur]
        overlaps = self.tden_overlaps(
            mo_coeffs_ref, ci_coeffs_ref, mo_coeffs_cur, ci_coeffs_cur, ao_ovlp
        )
        return overlaps

    def tden_overlap_with_calculator(self, calc, ao_ovlp=None):
        mo_coeffs1 = self.mo_coeff_list[-1]
        ci_coeffs1 = self.ci_coeff_list[-1]
        mo_coeffs2 = calc.mo_coeff_list[-1]
        ci_coeffs2 = calc.ci_coeff_list[-1]
        overlaps = self.tden_overlaps(
            mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2, ao_ovlp=ao_ovlp
        )
        return overlaps

    def calculate_state_ntos(self, state_ci_coeffs, mos):
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

    def get_nto_overlaps(self, indices=None, ao_ovlp=None, org=False):
        ref, cur = self.get_indices(indices)

        if ao_ovlp is None:
            ao_ovlp = self.get_sao_from_mo_coeffs_and_dump(
                self.get_ref_mos(self.mo_coeff_list[ref], self.mo_coeff_list[cur])
            )

        ntos_1 = self.nto_list[ref]
        ntos_2 = self.nto_list[cur]
        if org:
            overlaps = self.nto_org_overlaps(
                ntos_1, ntos_2, ao_ovlp, nto_thresh=self.nto_thresh
            )
        else:
            overlaps = self.nto_overlaps(ntos_1, ntos_2, ao_ovlp)
        return overlaps

    def nto_overlaps(self, ntos_1, ntos_2, ao_ovlp):
        states1 = len(ntos_1)
        states2 = len(ntos_2)
        ovlps = np.zeros((states1, states2))
        for i in range(states1):
            n_i = ntos_1[i]
            l_i = n_i.lambdas[:, None]
            ntos_i = l_i * n_i.ntos
            for j in range(i, states2):
                n_j = ntos_2[j]
                l_j = n_j.lambdas[:, None]
                ntos_j = l_j * n_j.ntos
                ovlp = np.sum(np.abs(ntos_i.dot(ao_ovlp).dot(ntos_j.T)))
                ovlps[i, j] = ovlp
                ovlps[j, i] = ovlp
        return ovlps

    def nto_org_overlaps(self, ntos_1, ntos_2, ao_ovlp, nto_thresh=0.3):
        states_1 = len(ntos_1)
        states_2 = len(ntos_2)
        ovlps = np.zeros((states_1, states_2))

        for i in range(states_1):
            n_i = ntos_1[i]
            l_i = n_i.lambdas[:, None]
            ntos_i = n_i.ntos[(l_i >= nto_thresh).flatten()]
            l_i_big = l_i[l_i >= nto_thresh]
            for j in range(i, states_2):
                n_j = ntos_2[j]
                l_j = n_j.lambdas[:, None]
                ntos_j = n_j.ntos[(l_j >= nto_thresh).flatten()]
                ovlp = np.sum(
                    l_i_big[:, None] * np.abs(ntos_i.dot(ao_ovlp).dot(ntos_j.T))
                )
                ovlps[i, j] = ovlp
                ovlps[j, i] = ovlp
        return ovlps

    def prepare_overlap_data(self, path):
        """Implement calculator specific parsing of MO coefficients and CI
        coefficients here. Should return a filename pointing to TURBOMOLE
        like mos, a MO coefficient array and a CI coefficient array."""
        raise Exception("Implement me!")

    def dump_overlap_data(self):
        if self.h5_dump:
            h5_group = self.get_h5_group()

            h5_group.attrs["ovlp_type"] = self.ovlp_type
            h5_group.attrs["ovlp_with"] = self.ovlp_with
            h5_group.attrs["orient"] = self.orient
            h5_group.attrs["atoms"] = np.string_(self.atoms)

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
            "mo_coeffs": np.array(self.mo_coeff_list, dtype=float),
            "ci_coeffs": np.array(self.ci_coeff_list, dtype=float),
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
                handle.create_dataset(name=key, dtype=val.dtype, data=val)
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
            mo_coeffs = handle["mo_coeffs"][:]
            ci_coeffs = handle["ci_coeffs"][:]
            all_energies = handle["all_energies"][:]
            try:
                ref_roots = handle["ref_roots"][:]
                roots = handle["roots"][:]
                calculated_roots = handle["calculated_roots"][:]
                root_info = True
            except KeyError:
                print(f"Couldn't find root information in '{h5_fn}'.")

        calc.ovlp_type = ovlp_type
        calc.ovlp_with = ovlp_with
        calc.mo_coeff_list = list(mo_coeffs)
        calc.ci_coeff_list = list(ci_coeffs)
        calc.all_energies_list = list(all_energies)
        if root_info:
            calc.roots_list = list(roots)
            calc.calculated_roots = list(calculated_roots)
            try:
                calc.first_root = ref_roots[0]
                calc.root = calc.first_root
            except IndexError:
                calc.root = roots[0]

        if (ovlp_type == "wf") or set_wfow:
            calc.set_wfow(ci_coeffs[0])

        return calc

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

    def update_array_dims(self):
        """Assure consistent array shapes when the number of calculated
        roots changed."""
        raise Exception("This method is not functional right now!")

    def set_wfow(self, ci_coeffs):
        occ_mo_num, virt_mo_num = ci_coeffs[0].shape
        try:
            wfow_mem = self.pal * self.mem
        except AttributeError:
            wfow_mem = 8000
        self.wfow = WFOWrapper(
            occ_mo_num,
            virt_mo_num,
            calc_number=self.calc_number,
            wfow_mem=wfow_mem,
            ncore=self.ncore,
            conf_thresh=self.conf_thresh,
        )

    def store_overlap_data(self, atoms, coords, path=None, overlap_data=None):
        if self.atoms is None:
            self.atoms = atoms

        if overlap_data is None:
            overlap_data = self.prepare_overlap_data(path)
        mo_coeffs, X, Y, all_ens = overlap_data

        ao_ovlp = self.get_sao_from_mo_coeffs(mo_coeffs)
        mo_coeffs = self.renorm_mos(mo_coeffs, ao_ovlp)
        if self.XY == "X":
            ci_coeffs = X
        elif self.XY == "X+Y":
            ci_coeffs = X + Y
        elif self.XY == "X-Y":
            ci_coeffs = X - Y
        else:
            raise Exception(
                f"Invalid 'XY' value. Allowed values are: '{self.VALID_XY}'!"
            )

        # Norm (X+Y) to 1 for every state
        ci_norms = np.linalg.norm(ci_coeffs, axis=(1, 2))
        ci_coeffs /= ci_norms[:, None, None]

        ci_norms = np.linalg.norm(ci_coeffs, axis=(1, 2))
        self.log(f"CI-vector norms: {ci_norms}")

        # Don't create the object when we use a different ovlp method.
        if (self.ovlp_type == "wf") and (self.wfow is None):
            self.set_wfow(ci_coeffs)

        if self.first_root is None:
            self.first_root = self.root
            self.log(f"Set first root to {self.first_root}.")
        # Used for transition density overlaps
        self.mo_coeff_list.append(mo_coeffs)
        self.ci_coeff_list.append(ci_coeffs)
        self.coords_list.append(coords)
        self.calculated_roots.append(self.root)
        # We can't calculate any overlaps in the first cycle, so we can't
        # compute a new root value. So we store the same value as for
        # calculated_roots.
        if len(self.ci_coeff_list) < 2:
            self.roots_list.append(self.root)
        self.all_energies_list.append(all_ens)
        # Also store NTOs if requested
        if self.ovlp_type in ("nto", "nto_org"):
            self.set_ntos(mo_coeffs, ci_coeffs)

        # self.update_array_dims()

    def track_root(self, ovlp_type=None):
        """Check if a root flip occured occured compared to the previous cycle
        by calculating the overlap matrix wrt. a reference cycle."""

        if ovlp_type is None:
            ovlp_type = self.ovlp_type

        # Nothing to compare to if only one calculation was done yet.
        # Nonetheless, dump the first cycle to HDF5.
        if len(self.ci_coeff_list) < 2:
            self.dump_overlap_data()
            self.log(
                "Skipping overlap calculation in the first cycle "
                "as there is nothing to compare to."
            )
            return False

        ao_ovlp = None
        # We can only run a double molecule calculation if it is
        # implemented for the specific calculator.
        if self.double_mol and hasattr(self, "run_double_mol_calculation"):
            old, new = self.get_indices()
            two_coords = self.coords_list[old], self.coords_list[new]
            ao_ovlp = self.run_double_mol_calculation(self.atoms, *two_coords)
        elif (self.double_mol is False) and (self.ovlp_type == "wf"):
            ao_ovlp = self.get_sao_from_mo_coeffs(self.mo_coeff_list[-1])
            self.log("Creating S_AO by myself to avoid its creation in " "WFOverlap.")

        if ovlp_type == "wf":
            overlap_mats = self.get_wf_overlaps(ao_ovlp=ao_ovlp)
            overlaps = np.abs(overlap_mats[2])
            # overlaps = overlaps**2
        elif ovlp_type == "tden":
            overlaps = self.get_tden_overlaps(ao_ovlp=ao_ovlp)
        elif ovlp_type == "nto":
            overlaps = self.get_nto_overlaps(ao_ovlp=ao_ovlp)
        elif ovlp_type == "nto_org":
            overlaps = self.get_nto_overlaps(ao_ovlp=ao_ovlp, org=True)
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
        self.log(
            f"Reference is cycle {self.ref_cycle}, root {ref_root}. "
            f"Analyzing row {row_ind} of the overlap matrix."
        )

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
        ci_coeffs = self.ci_coeff_list[cycle]
        exc_str = get_mwfn_exc_str(energies, ci_coeffs)
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
