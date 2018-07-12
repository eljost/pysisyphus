from pysisyphus.calculators.Calculator import Calculator


class OverlapCalculator(Calculator):

    def __init__(self, *args, track=False, ovlp_type="tden", **kwargs, ):
        self.track = track
        self.ovlp_type = ovlp_type
        self.mo_coeff_list = list()
        self.ci_coeff_list = list()

        super().__init__(*args, **kwargs)

    def tden_overlaps(self, mo_coeffs1, mo_coeffs2, ci_coeffs1, ci_coeffs2,
                      S_AO=None):
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
        S_AO : ndarray, shape(AOs1, AOs2)
            Double molcule AO overlaps.
        """

        mo_coeffs1_inv = np.linalg.inv(mo_coeffs1)
        # AO overlaps
        if S_AO is None:
            S_AO = mo_coeffs1_inv.dot(mo_coeffs1_inv.T)
        # MO overlaps
        S_MO = mo_coeffs1.dot(S_AO).dot(mo_coeffs2.T)

        overlaps = np.array(
            [np.sum(state1.dot(S_MO) * state2)
             for state1, state2 in it.product(ci_coeffs1, ci_coeffs2)
        ])

        return overlaps

    def track_tden_root(S_AO=None):
        mo_coeffs1 = self.mo_coeffs_list[-1]
        ci_coeffs1 = self.ci_coeff_list[-1]
        mo_coeffs2 = self.mo_coeffs_list[-2]
        ci_coeffs2 = self.ci_coeff_list[-2]
        overlaps = self.tden_overlaps(mo_coeffs1, ci_coeffs1,
                                      mo_coeffs2, ci_coeffs2,
                                      S_AO)
        raise Exception("wohoo")

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

    def track_root(self, atoms, coords, double_mol=True,
                   ovlp_type=None):
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
        if double_mol and hasattr(self, "run_double_mol_calculation"):
            last_two_coords = self.wfow.last_two_coords
            ao_ovlp = self.run_double_mol_calculation(atoms, *last_two_coords)

        if ovlp_type == "wf":
            self.root = self.wfow.track(old_root=self.root, ao_ovlp=ao_ovlp)
        elif ovlp_type == "tden":
            self.track_tden_root(ao_ovlp)
        else:
            raise Exception("Invalid overlap specifier! Use one of "
                            "'tden'/'wf'!")

        if self.root != old_root:
            self.log(f"Found a root flip from {old_root} to {self.root}!")

        # True if a root flip occured
        return not (self.root == old_root)
