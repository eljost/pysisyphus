#!/usr/bin/env python3

# [1] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b01148
#     Plasser, 2016
# [2] https://doi.org/10.1002/jcc.25800
#     Garcia, Campetella, 2019

from collections import namedtuple
from pathlib import Path
import shutil
import tempfile

import h5py
import itertools as it
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.WFOWrapper import WFOWrapper
from pysisyphus.constants import AU2EV
from pysisyphus.wrapper.mwfn import make_cdd
from pysisyphus.wrapper.jmol import render_cdd_cube as render_cdd_cube_jmol


NTOs = namedtuple("NTOs", "ntos lambdas")


class OverlapCalculator(Calculator):
    OVLP_TYPE_VERBOSE = {
        "wf": "wavefunction overlap",
        "tden": "transition density matrix overlap",
        "nto": "natural transition orbital overlap",
        # As described in 10.1002/jcc.25800
        "nto_org": "original natural transition orbital overlap",
    }
    VALID_KEYS = [k for k in OVLP_TYPE_VERBOSE.keys()]  # lgtm [py/non-iterable-in-for-loop]

    def __init__(self, *args, track=False, ovlp_type="wf", double_mol=False,
                 ovlp_with="previous", adapt_args=(0.5, 0.3, 0.6),
                 use_ntos=4, cdds=None, orient="", dump_fn="overlap_data.h5",
                 ncore=0, conf_thresh=1e-4, dyn_roots=0,
                 **kwargs):
        super().__init__(*args, **kwargs)

        self.track = track
        self.ovlp_type = ovlp_type
        assert self.ovlp_type in self.OVLP_TYPE_VERBOSE.keys(), \
		f"Valid overlap types are {self.VALID_KEYS}"
        self.double_mol = double_mol
        assert ovlp_with in ("previous", "first", "adapt")
        self.ovlp_with = ovlp_with
        self.adapt_args = np.abs(adapt_args, dtype=float)
        self.adpt_thresh, self.adpt_min, self.adpt_max = self.adapt_args
        self.use_ntos = use_ntos
        self.cdds = cdds
        self.orient = orient
        self.dump_fn = self.out_dir / dump_fn
        self.ncore = int(ncore)
        self.conf_thresh = float(conf_thresh)
        self.dyn_roots = int(dyn_roots)
        if self.dyn_roots != 0:
            self.dyn_roots = 0
            self.log("dyn_roots = 0 is hardcoded right now")

        assert self.ncore >= 0, "ncore must be a >= 0!"
        if self.cdds:
            assert self.cdds in "calc render".split()

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
        self.root_flips = [False, ]
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
            self.log("Tracking excited states with "
                    f"{self.OVLP_TYPE_VERBOSE[ovlp_type]}s "
                    f"between the current and the {self.ovlp_with} geometry."
            )
            if self.ovlp_with == "adapt":
                self.log(f"Adapt args: {self.adapt_args}")

    @property
    def roots_number(self):
        return self.root + self.dyn_roots

    def blowup_ci_coeffs(self, ci_coeffs):
        states, occ, virt = ci_coeffs.shape
        full = np.zeros((states, occ, occ+virt))
        full[:,:,occ:] = ci_coeffs
        return full

    def get_indices(self, indices=None):
        """
        A new root is determined by selecting the overlap matrix row
        corresponding to the old root and checking for the new root
        with the highest overlap (at the new geometry).

        The overlap matrix is usually formed by a double loop like:

        overlap_matrix = np.empty((old_states, new_states))
        for i, old_state in enumerate(old_states):
            for j, new_state in enumerate(new_states):
                overlap_matrix[i, j] = make_overlap(old_state, new_state)

        So the old states run along the rows. Thats why the old_state index
        comes first in the 'indices' tuple.
        """

        if indices is None:
            # By default we compare a reference cycle with the current (last)
            # cycle, so the second index is -1.
            old, new = self.ref_cycle, -1
        else:
            assert len(indices) == 2
            old, new = [int(i) for i in indices]
        return (old, new)

    def get_wf_overlaps(self, indices=None, ao_ovlp=None):
        old, new = self.get_indices(indices)
        old_cycle = (self.mo_coeff_list[old], self.ci_coeff_list[old])
        new_cycle = (self.mo_coeff_list[new], self.ci_coeff_list[new])
        return self.wfow.wf_overlap(old_cycle, new_cycle, ao_ovlp)

    def tden_overlaps(self, mo_coeffs1, ci_coeffs1, mo_coeffs2, ci_coeffs2,
                      ao_ovlp=None):
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
        states, occ, virt = ci_coeffs1.shape
        ci_full1 = self.blowup_ci_coeffs(ci_coeffs1)
        ci_full2 = self.blowup_ci_coeffs(ci_coeffs2)

        # AO overlaps
        if ao_ovlp is None:
            mo_coeffs2_inv = np.linalg.inv(mo_coeffs2)
            ao_ovlp = mo_coeffs2_inv.dot(mo_coeffs2_inv.T)
            ao_ovlp_fn = self.make_fn("ao_ovlp_rec")
            np.savetxt(ao_ovlp_fn, ao_ovlp)
        # MO overlaps
        S_MO = mo_coeffs1.dot(ao_ovlp).dot(mo_coeffs2.T)
        S_MO_occ = S_MO[:occ, :occ]

        overlaps = [np.sum(S_MO_occ.dot(state1).dot(S_MO) * state2)
                    for state1, state2 in it.product(ci_full1, ci_full2)
        ]
        overlaps = np.array(overlaps).reshape(states, -1)

        return overlaps

    def get_tden_overlaps(self, indices=None, ao_ovlp=None):
        old, new = self.get_indices(indices)
        mo_coeffs1 = self.mo_coeff_list[old]
        ci_coeffs1 = self.ci_coeff_list[old]
        mo_coeffs2 = self.mo_coeff_list[new]
        ci_coeffs2 = self.ci_coeff_list[new]
        overlaps = self.tden_overlaps(mo_coeffs1, ci_coeffs1,
                                      mo_coeffs2, ci_coeffs2,
                                      ao_ovlp=ao_ovlp)
        return overlaps

    def tdens_overlap_with_calculator(self, calc, ao_ovlp=None):
        mo_coeffs1 = self.mo_coeff_list[-1]
        ci_coeffs1 = self.ci_coeff_list[-1]
        mo_coeffs2 = calc.mo_coeff_list[-1]
        ci_coeffs2 = calc.ci_coeff_list[-1]
        overlaps = self.tden_overlaps(mo_coeffs1, ci_coeffs1,
                                      mo_coeffs2, ci_coeffs2,
                                      ao_ovlp=ao_ovlp)
        return overlaps

    def calculate_state_ntos(self, state_ci_coeffs, mos):
        normed = state_ci_coeffs / np.linalg.norm(state_ci_coeffs)
        # u, s, vh = np.linalg.svd(state_ci_coeffs)
        u, s, vh = np.linalg.svd(normed)
        lambdas = s**2
        self.log("Normalized transition density vector to 1.")
        self.log(f"Sum(lambdas)={np.sum(lambdas):.4f}")
        lambdas_str = np.array2string(lambdas[:3], precision=4,
                                      suppress_small=True)
        self.log(f"First three lambdas: {lambdas_str}")

        occ_mo_num = state_ci_coeffs.shape[0]
        occ_mos = mos[:occ_mo_num]
        vir_mos = mos[occ_mo_num:]
        occ_ntos = occ_mos.T.dot(u)
        vir_ntos = vir_mos.T.dot(vh)
        return occ_ntos, vir_ntos, lambdas

    def get_nto_overlaps(self, indices=None, ao_ovlp=None, org=False):
        old, new = self.get_indices(indices)

        if ao_ovlp is None:
            mos_inv = np.linalg.inv(self.mo_coeff_list[new])
            ao_ovlp = mos_inv.dot(mos_inv.T)

        ntos_1 = self.nto_list[old]
        ntos_2 = self.nto_list[new]
        if org:
            overlaps = self.nto_org_overlaps(ntos_1, ntos_2, ao_ovlp)
        else:
            overlaps = self.nto_overlaps(ntos_1, ntos_2, ao_ovlp)
        return overlaps

    def nto_overlaps(self, ntos_1, ntos_2, ao_ovlp):
        states1 = len(ntos_1)
        states2 = len(ntos_2)
        ovlps = np.zeros((states1, states2))
        for i in range(states1):
            n_i = ntos_1[i]
            l_i = n_i.lambdas[:,None]
            ntos_i = l_i*n_i.ntos
            for j in range(i, states2):
                n_j = ntos_2[j]
                l_j = n_j.lambdas[:,None]
                ntos_j = l_j*n_j.ntos
                ovlp = np.sum(np.abs(ntos_i.dot(ao_ovlp).dot(ntos_j.T)))
                ovlps[i, j] = ovlp
                ovlps[j, i] = ovlp
        return ovlps

    def nto_org_overlaps(self, ntos_1, ntos_2, ao_ovlp, nto_thresh=.3):
        states_1 = len(ntos_1)
        states_2 = len(ntos_2)
        ovlps = np.zeros((states_1, states_2))

        for i in range(states_1):
            n_i = ntos_1[i]
            l_i = n_i.lambdas[:,None]
            ntos_i = n_i.ntos[(l_i >= nto_thresh).flatten()]
            l_i_big = l_i[l_i >= nto_thresh]
            for j in range(i, states_2):
                n_j = ntos_2[j]
                l_j = n_j.lambdas[:,None]
                ntos_j = n_j.ntos[(l_j >= nto_thresh).flatten()]
                ovlp = np.sum(l_i_big[:,None] * np.abs(ntos_i.dot(ao_ovlp).dot(ntos_j.T)))
                ovlps[i, j] = ovlp
                ovlps[j, i] = ovlp
        return ovlps

    def prepare_overlap_data(self):
        """Implement calculator specific parsing of MO coefficients and CI
        coefficients here. Should return a filename pointing to TURBOMOLE
        like mos, a MO coefficient array and a CI coefficient array."""
        raise Exception("Implement me!")

    def dump_overlap_data(self):
        data_dict = {
            "mo_coeffs": np.array(self.mo_coeff_list, dtype=float),
            "ci_coeffs": np.array(self.ci_coeff_list, dtype=float),
            "coords": np.array(self.coords_list, dtype=float),
            "all_energies": np.array(self.all_energies_list, dtype=float),
            "orient": np.array(self.orient, dtype="S"),
            "atoms": np.array(self.atoms, dtype="S")
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
            handle.create_dataset(name="ovlp_type", data=np.string_(self.ovlp_type))
            handle.create_dataset(name="ovlp_with", data=np.string_(self.ovlp_with))

    @staticmethod
    def from_overlap_data(h5_fn):
        calc_kwargs = {
            "track": True,
            # "ovlp_with": ovlp_with,
            # "ovlp_type": ovlp_type,
        }
        calc = OverlapCalculator(**calc_kwargs)

        root_info = False
        with h5py.File(h5_fn) as handle:
            ovlp_with = handle["ovlp_with"][()].decode()
            ovlp_type = handle["ovlp_type"][()].decode()
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

        calc.mo_coeff_list = list(mo_coeffs)
        calc.ci_coeff_list = list(ci_coeffs)
        calc.all_energies_list = list(all_energies)
        if root_info:
            calc.roots_list = list(roots)
            calc.calculated_roots = list(calculated_roots)
            calc.first_root = ref_roots[0]
            calc.root = calc.first_root

        return calc

    def set_ntos(self, mo_coeffs, ci_coeffs):
        roots = ci_coeffs.shape[0]
        ntos_for_cycle = list()
        for root in range(roots):
            sn_ci_coeffs = ci_coeffs[root]
            self.log("Calculating NTOs for root {root+1}")
            occ_ntos, vir_ntos, lambdas = self.calculate_state_ntos(
                                                    sn_ci_coeffs,
                                                    mo_coeffs,
            )
            ovlp_occ_ntos = occ_ntos.T[:self.use_ntos]
            ovlp_vir_ntos = vir_ntos.T[:self.use_ntos]
            ovlp_lambdas = lambdas[:self.use_ntos]
            ovlp_lambdas = np.concatenate((ovlp_lambdas, ovlp_lambdas))
            ovlp_ntos = np.concatenate((ovlp_occ_ntos, ovlp_vir_ntos), axis=0)
            ntos = NTOs(ntos=ovlp_ntos, lambdas=ovlp_lambdas)
            ntos_for_cycle.append(ntos)
        self.nto_list.append(ntos_for_cycle)

    def update_array_dims(self):
        """Assure consistent array shapes when the number of calculated
        roots changed."""
        raise Exception("This method is not functional right now!")

        # # We only have to update the first dimension
        # to_update = (
            # "all_energies_list",
            # "overlap_matrices",
            # "ci_coeff_list",
        # )
        # for name in to_update:
            # arr_list = getattr(self, name)
            # first_dims = [arr.shape[0] for arr in arr_list]
            # max_ = max(first_dims)
            # print(f"name: {name} first_dims: {first_dims}")

        # try:
            # # With dyn_root > 0 the size of all_ens may be different from cycle to
            # # cycle. Here we assure that the stored array will always have the same
            # # size and shape as the first array in the list.
            # all_ens_zero = np.full_like(self.all_energies_list[0], np.nan, dtype=float)
            # all_ens_zero[:all_ens.size] = all_ens
            # all_ens = all_ens_zero
        # except IndexError:
            # pass

    def store_overlap_data(self, atoms, coords):
        if self.atoms is None:
            self.atoms = atoms
        mo_coeffs, ci_coeffs, all_ens = self.prepare_overlap_data()

        # Don't create the object when we use a different ovlp method.
        if (self.ovlp_type == "wf") and (self.wfow is None):
            occ_mo_num, virt_mo_num = ci_coeffs[0].shape
            try:
                wfow_mem = self.pal * self.mem
            except AttributeError:
                wfow_mem = 8000
            self.wfow = WFOWrapper(occ_mo_num, virt_mo_num,
                                   calc_number=self.calc_number, wfow_mem=wfow_mem,
                                   ncore=self.ncore, conf_thresh=self.conf_thresh,
            )

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
        # Nothing to compare to if only one calculation was done yet
        if len(self.ci_coeff_list) < 2:
            self.log("Skipping overlap calculation in the first cycle "
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
            mos_inv = np.linalg.inv(self.mo_coeff_list[-1])
            ao_ovlp = mos_inv.dot(mos_inv.T)
            self.log("Creating S_AO by myself to avoid its creation in "
                     "WFOverlap.")

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
            raise Exception("Invalid overlap type key! Use one of "
			    + ", ".join(self.VALID_KEYS))
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
        self.log(f"Reference is cycle {self.ref_cycle}, root {ref_root}. "
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
            self.log(f"Root flip! New root is {self.root}. Root at previous "
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
            self.log(f"Two highest overlaps: {sec_highest:.2%}, {highest:.2%}, "
                     f"ratio={ratio:.4f}")
            above_thresh = highest >= self.adpt_thresh
            self.log(f"Highest overlap is above threshold? (>= {self.adpt_thresh:.4f}): "
                     f"{above_thresh}")
            valid_ratio = self.adpt_min < ratio < self.adpt_max
            self.log(f"Ratio is valid? (between {self.adpt_min:.4f} and "
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
            self.calc_cdd_cube(self.root)
            if self.cdds == "render":
                self.render_cdd_cube()

        self.dump_overlap_data()

        # True if a root flip occured
        return root_flip

    def get_mwfn_exc_str(self, cycle, thresh=1e-3):
        energies = self.all_energies_list[cycle]
        # Cycle, states, occ, virt
        ci_coeffs = self.ci_coeff_list[cycle]
        exc_energies = (energies[1:] - energies[0]) * AU2EV
        above_thresh = np.abs(ci_coeffs) > thresh
        _, occ_mos, virt_mos = ci_coeffs.shape

        exc_str = ""
        mult = 1
        self.log(f"Using dummy multiplicity={mult} in get_mwfn_exc_str")
        for root_, (root_ci_coeffs, exc_en) in enumerate(zip(ci_coeffs, exc_energies), 1):
            exc_str += f"Excited State {root_} {mult} {exc_en:.4f}\n"
            for (occ, virt), coeff in np.ndenumerate(root_ci_coeffs):
                if abs(coeff) < thresh:
                    continue
                occ_mo = occ+1
                virt_mo = occ_mos + 1 + virt
                exc_str += f"{occ_mo:>8d} -> {virt_mo}       {coeff: .5f}\n"
            exc_str += "\n"
        return exc_str

    def calc_cdd_cube(self, root, cycle=-1):
        if (not hasattr(self, "mwfn_wf")):
            self.log("self.mwfn_wf is not set! Skipping CDD cube generation!")
        if cycle != -1:
            self.log("'cycle' argument to make_cdd_cube is currently ignored!")
        exc_str = self.get_mwfn_exc_str(cycle)
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
