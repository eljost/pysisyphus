#!/usr/bin/env python3

# [1] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b01148
#     Plasser, 2016

from collections import namedtuple
from pathlib import Path
import shutil
import tempfile

import h5py
import itertools as it
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2EV
from pysisyphus.wrapper.mwfn import make_cdd


NTOs = namedtuple("NTOs", "ntos lambdas")


class OverlapCalculator(Calculator):
    OVLP_TYPE_VERBOSE = {
        "wf": "wavefunction overlap",
        "tden": "transition density matrix overlap",
        "nto": "natural transition orbital overlap",
        # As described in 10.1002/jcc.25800
        "nto_org": "original natural transition orbital overlap",
    }
    VALID_KEYS = [k for k in OVLP_TYPE_VERBOSE.keys()]

    def __init__(self, *args, track=False, ovlp_type="wf", double_mol=False,
                 ovlp_with="previous", adapt_args=(0.5, 0.3, 0.6),
                 use_ntos=4, make_cubes=False, **kwargs):
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
        self.make_cubes = make_cubes

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
        self.dump_fn = "overlap_data.h5"

        if track:
            self.log("Tracking excited states with "
                    f"{self.OVLP_TYPE_VERBOSE[ovlp_type]}s "
                    f"between the current and the {self.ovlp_with} geometry."
            )
            if self.ovlp_with == "adapt":
                self.log(f"Adapt args: {self.adapt_args}")

    def blowup_ci_coeffs(self, ci_coeffs):
        states, occ, virt = ci_coeffs.shape
        full = np.zeros((states, occ, occ+virt))
        full[:,:,occ:] = ci_coeffs
        return full

    def get_indices(self):
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

        # ref_inds = {
            # "first": 0,
            # # Data of the current cycle is saved before this method is called,
            # # so the current cycle resides at -1 and the previous cycle is at -2.
            # "previous": -2,
            # "adapt": self.ref_cycle,
        # }
        # ref_ind = ref_inds[self.ovlp_with]
        # ref_ind =

        # We always compare against the current cycle, so the second index
        # is always -1.
        return (self.ref_cycle, -1)

    def get_wfow_overlaps(self, ao_ovlp=None):
        old, new = self.get_indices()
        return self.wfow.overlaps(old, new, ao_ovlp)

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
            mo_coeffs1_inv = np.linalg.inv(mo_coeffs1)
            ao_ovlp = mo_coeffs1_inv.dot(mo_coeffs1_inv.T)
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

    def get_tden_overlaps(self, ao_ovlp=None):
        old, new = self.get_indices()
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

    def get_nto_overlaps(self, ao_ovlp=None, org=False):
        if ao_ovlp is None:
            mos_inv = np.linalg.inv(self.mo_coeff_list[-1])
            ao_ovlp = mos_inv.dot(mos_inv.T)

        old, new = self.get_indices()
        ntos_1 = self.nto_list[old]
        ntos_2 = self.nto_list[new]
        if org:
            overlaps = self.nto_org_overlaps(ntos_1, ntos_2, ao_ovlp)
        else:
            overlaps = self.nto_overlaps(ntos_1, ntos_2, ao_ovlp)
        return overlaps

    def nto_overlaps(self, ntos_1, ntos_2, ao_ovlp):
        assert len(ntos_1) == len(ntos_2), "For now only the same number of "\
                                           "states is required."
        states = len(ntos_1)
        ovlps = np.zeros((states, states))
        for i in range(states):
            n_i = ntos_1[i]
            l_i = n_i.lambdas[:,None]
            ntos_i = l_i*n_i.ntos
            for j in range(i, states):
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
            "calculated_roots": np.array(self.calculated_roots, dtype=int),
            "roots": np.array(self.roots_list, dtype=int),
            "all_energies": np.array(self.all_energies_list, dtype=float),
            "root_flips": np.array(self.root_flips, dtype=bool),
            "overlap_matrices": np.array(self.overlap_matrices, dtype=float),
            "row_inds": np.array(self.row_inds, dtype=int),
            "ref_cycles": np.array(self.ref_cycles, dtype=int),
        }

        with h5py.File(self.dump_fn, "w") as handle:
            for key, val in data_dict.items():
                handle.create_dataset(name=key, dtype=val.dtype, data=val)
            handle.create_dataset(name="ovlp_type", data=np.string_(self.ovlp_type))
            handle.create_dataset(name="ovlp_with", data=np.string_(self.ovlp_with))

    def store_overlap_data(self, atoms, coords):
        mos_fn, mo_coeffs, ci_coeffs, all_ens = self.prepare_overlap_data()
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
        # Used for WFOverlap
        self.wfow.store_iteration(atoms, coords, mos_fn, ci_coeffs)
        # Also store NTOs if requested
        self.all_energies_list.append(all_ens)
        if self.ovlp_type in ("nto", "nto_org"):
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

    def track_root(self, atoms, coords, ovlp_type=None):
        """Store the information of the current iteration and if possible
        calculate the overlap with a previous/the first iteration."""
        self.store_overlap_data(atoms, coords)
        # Calculated root and self.root are still the same right now.
        assert self.calculated_roots == self.root
        if self.make_cubes and len(self.ci_coeff_list >= 1):
            self.make_cdd_cube(self.root, cycle=-1)

        old_root = self.root
        if not ovlp_type:
            ovlp_type = self.ovlp_type
        # Nothing to compare to if only one calculation was done yet
        if len(self.ci_coeff_list) < 2:
            self.log("Skipping overlap calculation in the first cycle "
                     "as there is nothing to compare to."
            )
            return False

        ao_ovlp = None
        # We can only run a double molecule calculation if it is
        # implemented for the specific calculator, so we have to check it.
        if self.double_mol and hasattr(self, "run_double_mol_calculation"):
            old, new = self.get_indices()
            two_coords = self.coords_list[old], self.coords_list[new]
            ao_ovlp = self.run_double_mol_calculation(atoms, *two_coords)
        elif (self.double_mol is False) and (self.ovlp_type == "wf"):
            mos_inv = np.linalg.inv(self.mo_coeff_list[-1])
            ao_ovlp = mos_inv.dot(mos_inv.T)
            self.log("Creating S_AO by myself to avoid its creation in "
                     "WFOverlap.")

        if ovlp_type == "wf":
            overlap_mats = self.get_wfow_overlaps(ao_ovlp)
            overlaps = np.abs(overlap_mats[2])
            # overlaps = overlaps**2
        elif ovlp_type == "tden":
            overlaps = self.get_tden_overlaps(ao_ovlp)
        elif ovlp_type == "nto":
            overlaps = self.get_nto_overlaps(ao_ovlp)
        elif ovlp_type == "nto_org":
            overlaps = self.get_nto_overlaps(ao_ovlp, org=True)
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
        self.dump_overlap_data()

        # True if a root flip occured
        return root_flip

    def get_mwfn_exc_str(self, cycle, thresh=1e-3):
        energies = self.all_energies_list[cycle]
        # Cycle, states, occ, virt
        ci_coeffs = self.ci_coeff_list[cycle]
        exc_energies = (energies[1:] - energies[0]) * AU2EV
        above_thresh = np.abs(ci_coeffs) > thresh
        occ_mos, virt_mos = ci_coeffs.shape[:2]

        exc_str = ""
        mult = 1
        self.log("Assuming mult={mult} in get_mwfn_exc_str")
        for root_, (ci_coeffs, exc_en) in enumerate(zip(ci_coeffs, exc_energies), 1):
            exc_str += f"Excited State {root_} {mult} {exc_en:.4f}\n"
            for (occ, virt), coeff in np.ndenumerate(ci_coeffs):
                if abs(coeff) < thresh:
                    continue
                occ_mo = occ+1
                virt_mo = occ_mos + 1 + virt
                exc_str += f"{occ_mo:>8d} -> {virt_mo}       {coeff: .5f}\n"
            exc_str += "\n"
        return exc_str

    def make_cdd_cube(self, root, cycle=-1):
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
            shutil.copy(cube, self.make_fn(cube_fn))
