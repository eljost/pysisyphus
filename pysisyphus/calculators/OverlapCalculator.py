from pysisyphus.calculators.Calculator import Calculator

class OverlapCalculator(Calculator):

    def track_root(self, atoms, coords, double_mol=True):
        """Store the information of the current iteration and if possible
        calculate the overlap with the previous iteration."""
        self.store_wfo_data(atoms, coords)
        old_root = self.root
        # Nothing to compare to if only one calculation was done yet
        if self.calc_counter <= 1:
            return False

        ao_ovlp = None
        if double_mol and hasattr(self, "run_double_mol_calculation"):
            last_two_coords = self.wfow.last_two_coords
            ao_ovlp = self.run_double_mol_calculation(atoms, *last_two_coords)

        self.root = self.wfow.track(old_root=self.root, ao_ovlp=ao_ovlp)
        if self.root != old_root:
            self.log(f"Found a root flip from {old_root} to {self.root}!")

        # True if a root flip occured
        return not (self.root == old_root)

    def set_mo_coeffs(self):
        raise Exception("Implement me!")

    def set_ci_coeffs(self):
        raise Exception("Implement me!")
