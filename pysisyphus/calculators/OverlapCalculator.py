import itertools as it
import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class OverlapCalculator(Calculator):
    ovlp_type_verbose = {
        "wf": "wavefunction overlap",
        "tden": "transition densisty matrix overlap",
    }


    def __init__(self, *args, track=False, ovlp_type="wf", double_mol=False,
                 **kwargs, ):
        self.track = track
        self.ovlp_type = ovlp_type
        self.double_mol = double_mol
        self.mo_coeff_list = list()
        self.ci_coeff_list = list()

        super().__init__(*args, **kwargs)

        if track:
            self.log("Tracking excited states with "
                    f"{self.ovlp_type_verbose[ovlp_type]}s")

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

    def last_two_tdens_overlap(self, ao_ovlp=None):
        mo_coeffs1 = self.mo_coeff_list[-1]
        ci_coeffs1 = self.ci_coeff_list[-1]
        mo_coeffs2 = self.mo_coeff_list[-2]
        ci_coeffs2 = self.ci_coeff_list[-2]
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
        # Used for transition density overlaps
        self.mo_coeff_list.append(mo_coeffs)
        self.ci_coeff_list.append(ci_coeffs)
        # Used for WFOverlap
        self.wfow.store_iteration(atoms, coords, mos_fn, ci_coeffs)


    def track_root(self, atoms, coords, ovlp_type=None):
        """Store the information of the current iteration and if possible
        calculate the overlap with the previous iteration."""
        self.store_overlap_data(atoms, coords)
        old_root = self.root
        if not ovlp_type:
            ovlp_type = self.ovlp_type
        # Nothing to compare to if only one calculation was done yet
        if self.calc_counter <= 1:
            return False

        ao_ovlp = None
        if self.double_mol and hasattr(self, "run_double_mol_calculation"):
            last_two_coords = self.wfow.last_two_coords
            ao_ovlp = self.run_double_mol_calculation(atoms, *last_two_coords)

        if ovlp_type == "wf":
            overlap_mats = self.wfow.overlaps(ao_ovlp)
            overlaps = overlap_mats[0]
            # overlaps = overlaps**2
        elif ovlp_type == "tden":
            overlaps = self.last_two_tdens_overlap(ao_ovlp)
            overlaps = np.abs(overlaps)
        else:
            raise Exception("Invalid overlap specifier! Use one of "
                            "'tden'/'wf'!")

        old_root = self.root
        self.log(f"Previous root is {old_root}.")
        old_root_col = overlaps[old_root-1]
        new_root = old_root_col.argmax()
        max_overlap = old_root_col[new_root]
        self.root = new_root + 1
        old_root_col_str = ", ".join(
            [f"{i}: {ov:.2%}" for i, ov in enumerate(old_root_col)]
        )
        self.log(f"Overlaps: {old_root_col_str}")
        root_flip = self.root != old_root
        if not root_flip:
            msg = f"New root is {self.root}, keeping previous root. Overlap is " \
                  f"{max_overlap:.2%}."
        else:
            msg = f"Root flip! New root {self.root} has {max_overlap:.2%} " \
                  f"overlap with previous root {old_root}."
        self.log(msg)

        # True if a root flip occured
        return root_flip
