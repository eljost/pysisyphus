import itertools as it
import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class OverlapCalculator(Calculator):
    ovlp_type_verbose = {
        "wf": "wavefunction overlap",
        "tden": "transition densisty matrix overlap",
    }

    def __init__(self, *args, track=False, ovlp_type="wf", double_mol=False,
                 ovlp_with="previous", **kwargs):
        self.track = track
        self.ovlp_type = ovlp_type
        self.double_mol = double_mol
        assert ovlp_with in ("previous", "first")
        self.ovlp_with = ovlp_with

        self.mo_coeff_list = list()
        self.ci_coeff_list = list()
        self.coords_list = list()
        self.first_root = None

        super().__init__(*args, **kwargs)

        if track:
            self.log("Tracking excited states with "
                    f"{self.ovlp_type_verbose[ovlp_type]}s using overlaps "
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

    def prepare_overlap_data(self):
        """Implement calculator specific parsing of MO coefficients and CI
        coefficients here. Should return a filename pointing to TURBOMOLE
        like mos, a MO coefficient array and a CI coefficient array."""
        raise Exception("Implement me!")

    def store_overlap_data(self, atoms, coords):
        mos_fn, mo_coeffs, ci_coeffs = self.prepare_overlap_data()
        if self.first_root is None:
            self.first_root = self.root
            self.log(f"Set first root to root {self.first_root}.")
        # Used for transition density overlaps
        self.mo_coeff_list.append(mo_coeffs)
        self.ci_coeff_list.append(ci_coeffs)
        self.coords_list.append(coords)
        # Used for WFOverlap
        self.wfow.store_iteration(atoms, coords, mos_fn, ci_coeffs)

    def track_root(self, atoms, coords, ovlp_type=None):
        """Store the information of the current iteration and if possible
        calculate the overlap with a previous/the first iteration."""
        self.store_overlap_data(atoms, coords)
        old_root = self.root
        if not ovlp_type:
            ovlp_type = self.ovlp_type
        # Nothing to compare to if only one calculation was done yet
        if self.calc_counter <= 1:
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

        if ovlp_type == "wf":
            overlap_mats = self.wfow.overlaps(ao_ovlp)
            overlaps = overlap_mats[0]
            # overlaps = overlaps**2
        elif ovlp_type == "tden":
            overlaps = self.get_tden_overlaps(ao_ovlp)
            overlaps = np.abs(overlaps)
        else:
            raise Exception("Invalid overlap specifier! Use one of "
                            "'tden'/'wf'!")

        prev_root = self.root
        self.log(f"Previous root is {prev_root}.")
        if self.ovlp_with == "first":
            prev_root_col = overlaps[self.first_root-1]
        elif self.ovlp_with == "previous":
            prev_root_col = overlaps[prev_root-1]
        new_root = prev_root_col.argmax()
        max_overlap = prev_root_col[new_root]
        self.root = new_root + 1
        prev_root_col_str = ", ".join(
            [f"{i}: {ov:.2%}" for i, ov in enumerate(prev_root_col)]
        )
        self.log(f"Overlaps: {prev_root_col_str}")
        root_flip = self.root != prev_root
        self.log(f"Highest overlap is {max_overlap:.2%}.")
        if not root_flip:
            self.log(f"Keeping current root {self.root}.")
        else:
            self.log(f"Root flip! New root is {self.root}. Root at previous "
                     f"step was {prev_root}."
            )

        # True if a root flip occured
        return root_flip
