#!/usr/bin/env python3

from collections import namedtuple

import h5py
import itertools as it
import numpy as np

from pysisyphus.calculators.Calculator import Calculator


NTOs = namedtuple("NTOs", "ntos lambdas")


class OverlapCalculator(Calculator):
    ovlp_type_verbose = {
        "wf": "wavefunction overlap",
        "tden": "transition density matrix overlap",
        "nto": "natural transition orbital overlap",
        # As described in 10.1002/jcc.25800
        "nto_org": "original natural transition orbital overlap",
    }

    def __init__(self, *args, track=False, ovlp_type="wf", double_mol=False,
                 ovlp_with="previous", use_ntos=4, **kwargs):
        self.track = track
        self.ovlp_type = ovlp_type
        assert self.ovlp_type in ("tden", "wf", "nto", "nto_org")
        self.double_mol = double_mol
        assert ovlp_with in ("previous", "first")
        self.ovlp_with = ovlp_with
        self.use_ntos = use_ntos

        self.mo_coeff_list = list()
        self.ci_coeff_list = list()
        self.nto_list = list()
        self.coords_list = list()
        self.roots_list = list()
        self.all_energies_list = list()
        self.root_flips_list = [False, ]
        self.first_root = None

        self.dump_fn = "overlap_data.h5"

        super().__init__(*args, **kwargs)

        if track:
            self.log("Tracking excited states with "
                    f"{self.ovlp_type_verbose[ovlp_type]}s "
                    f"between the current and the {self.ovlp_with} geometry."
            )

    def blowup_ci_coeffs(self, ci_coeffs):
        states, occ, virt = ci_coeffs.shape
        full = np.zeros((states, occ, occ+virt))
        full[:,:,occ:] = ci_coeffs
        return full

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
        # Overlap with previous cycle is the default
        indices = (-1, -2)
        if self.ovlp_with == "first":
            indices = (-1, 0)
        cur, prev = indices
        mo_coeffs1 = self.mo_coeff_list[cur]
        ci_coeffs1 = self.ci_coeff_list[cur]
        mo_coeffs2 = self.mo_coeff_list[prev]
        ci_coeffs2 = self.ci_coeff_list[prev]
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

        indices = (-1, -2)
        if self.ovlp_with == "first":
            indices = (-1, 0)
        ind_1, ind_2 = indices
        ntos_1 = self.nto_list[ind_1]
        ntos_2 = self.nto_list[ind_2]
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
            "roots": np.array(self.roots_list, dtype=int),
            "all_energies": np.array(self.all_energies_list, dtype=float),
            "root_flips": np.array(self.root_flips_list, dtype=bool),
        }
        # if self.nto_list:
            # data_dict["ntos"] = self.nto_list

        # with open(self.dump_fn, "w") as handle:
            # yaml.dump(data_dict, handle)
        # self.log(f"Saved OverlapCalculator data to '{self.dump_fn}'")
        with h5py.File(self.dump_fn, "w") as handle:
            for key, val in data_dict.items():
                handle.create_dataset(name=key, dtype=val.dtype, data=val)

    def store_overlap_data(self, atoms, coords):
        mos_fn, mo_coeffs, ci_coeffs, all_ens = self.prepare_overlap_data()
        if self.first_root is None:
            self.first_root = self.root
            self.log(f"Set first root to {self.first_root}.")
        # Used for transition density overlaps
        self.mo_coeff_list.append(mo_coeffs)
        self.ci_coeff_list.append(ci_coeffs)
        self.coords_list.append(coords)
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
        old_root = self.root
        if not ovlp_type:
            ovlp_type = self.ovlp_type
        # Nothing to compare to if only one calculation was done yet
        if len(self.ci_coeff_list) < 2:
            return False

        ao_ovlp = None
        # We can only run a double molecule calculation if it is
        # implemented for the specific calculator, so we have to check it.
        if self.double_mol and hasattr(self, "run_double_mol_calculation"):
            # Overlap with previous cycle is the default
            indices = (-1, -2)
            if self.ovlp_with == "first":
                indices = (-1, 0)
            cur, prev = indices
            two_coords = self.coords_list[cur], self.coords_list[prev]
            ao_ovlp = self.run_double_mol_calculation(atoms, *two_coords)
        elif (self.double_mol is False) and (self.ovlp_type == "wf"):
            mos_inv = np.linalg.inv(self.mo_coeff_list[-1])
            ao_ovlp = mos_inv.dot(mos_inv.T)
            self.log("Creating S_AO by myself to avoid its creation in "
                     "WFOverlap.")

        if ovlp_type == "wf":
            overlap_mats = self.wfow.overlaps(ao_ovlp)
            overlaps = overlap_mats[0]
            # overlaps = overlaps**2
        elif ovlp_type == "tden":
            overlaps = self.get_tden_overlaps(ao_ovlp)
            overlaps = np.abs(overlaps)
        elif ovlp_type == "nto":
            overlaps = self.get_nto_overlaps(ao_ovlp)
            overlaps = np.abs(overlaps)
        elif ovlp_type == "nto_org":
            overlaps = self.get_nto_overlaps(ao_ovlp, org=True)
            overlaps = np.abs(overlaps)
        else:
            raise Exception("Invalid overlap specifier! Use one of "
                            "'tden'/'wf'/'nto'!")

        prev_root = self.root
        self.log(f"Previous root is {prev_root}.")
        if self.ovlp_with == "first":
            row_ind = self.first_root-1
        elif self.ovlp_with == "previous":
            row_ind = prev_root - 1
        # With WFOverlaps the ground state is also present and the overlap
        # matrix has an additional row.
        if self.ovlp_type == "wf":
            row_ind += 1

        prev_root_row = overlaps[row_ind]
        new_root = prev_root_row.argmax()
        max_overlap = prev_root_row[new_root]
        if self.ovlp_type == "wf":
            new_root -= 1
        self.root = new_root + 1
        prev_root_row_str = ", ".join(
            [f"{i}: {ov:.2%}" for i, ov in enumerate(prev_root_row)]
        )
        self.log(f"Overlaps: {prev_root_row_str}")
        root_flip = self.root != prev_root
        self.log(f"Highest overlap is {max_overlap:.2%}.")
        if not root_flip:
            self.log(f"Keeping current root {self.root}.")
        else:
            self.log(f"Root flip! New root is {self.root}. Root at previous "
                     f"step was {prev_root}."
            )

        self.root_flips_list.append(root_flip)
        self.dump_overlap_data()

        # True if a root flip occured
        return root_flip
