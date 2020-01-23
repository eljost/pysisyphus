#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.guess_hessians import ts_hessian


class TSHessianOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    def __init__(self, geometry, root=0, hessian_ref=None, prim_coord=None,
                 rx_coords=None, hessian_init="calc", hessian_update="bofill",
                 max_micro_cycles=50, trust_radius=0.3, **kwargs):

        assert hessian_update == "bofill", \
            "Bofill update is recommended in a TS-optimization."

        super().__init__(geometry, hessian_init=hessian_init,
                         hessian_update=hessian_update, trust_radius=trust_radius,
                         **kwargs)

        self.root = int(root)
        self.hessian_ref = hessian_ref
        try:
            self.hessian_ref = np.loadtxt(self.hessian_ref)
            _ = self.geometry.coords.size
            expected_shape = (_, _)
            shape = self.hessian_ref.shape
            assert shape == expected_shape, \
                f"Shape of reference hessian {shape} doesn't match the expected "  \
                f"shape {expected_shape} of the hessian for the current coordinates!"
        except OSError as err:
            self.log(f"Tried to load reference hessian from '{self.hessian_ref}' "
                      "but the file could not be found.")
        except ValueError as err:
            self.log(f"No reference hessian provided.")
        self.prim_coord = prim_coord
        self.rx_coords = rx_coords

        self.ts_mode = None
        self.max_micro_cycles = max_micro_cycles
        self.prim_contrib_thresh = 0.05

        self.alpha0 = 1

    def prepare_opt(self):
        super().prepare_opt()

        # Assume a guess hessian when not calculated. This hessian has to be
        # modified according to the assumed reaction coordinates.
        if self.hessian_init != "calc":
            assert self.geometry.coord_type != "cart", \
                "Using a modified guess hessian for TS-optimizations is " \
                "only supported in redundand internal coordinates " \
                "(coord_type=redund)"
            prim_inds = [self.geometry.internal.get_index_of_prim_coord(rxc)
                         for rxc in self.rx_coords]
            missing_prim_inds = [self.rx_coords[i] for i, _ in enumerate(prim_inds)
                                 if _ is None]
            assert len(missing_prim_inds) == 0, \
                 "Some of the requested reaction coordinates are not defined: " \
                f"{missing_prim_inds}"
            self.H = ts_hessian(self.H, coord_inds=prim_inds)

        # Determiniation of initial mode either by using a provided
        # reference hessian, or by using a supplied root.

        eigvals, eigvecs = np.linalg.eigh(self.H)
        neg_inds = eigvals < -self.small_eigval_thresh

        self.log_negative_eigenvalues(eigvals, "Initial ")

        self.log("Determining initial TS mode to follow uphill.")
        # Select an initial TS-mode by highest overlap with eigenvectors from
        # reference hessian.
        if self.prim_coord:
            prim_ind = self.geometry.internal.get_index_of_prim_coord(self.prim_coord)
            if prim_ind is None:
                raise Exception(f"Primitive internal {self.prim_coord} is not defined!")
            # Select row of eigenvector-matrix that corresponds to this coordinate
            prim_row = eigvecs[prim_ind]
            big_contribs = np.abs(prim_row) > self.prim_contrib_thresh
            # Only consider negative eigenvalues
            big_contribs = np.bitwise_and(big_contribs, neg_inds)
            # Holds the indices of the modes to consider
            big_inds = np.arange(prim_row.size)[big_contribs]
            try:
                max_contrib_ind = big_inds[np.abs(prim_row[big_contribs]).argmax()]
            except ValueError as err:
                print( "No imaginary mode with significant contribution "
                      f"(>{self.prim_contrib_thresh:.3f}) of primitive internal "
                      f"{self.prim_coord} found!")
                raise err
            self.root = max_contrib_ind

            contrib_str = "\n".join(
                [f"\t{ind:03d}: {contrib:.4f}"
                 for ind, contrib in zip(big_inds, np.abs(prim_row)[big_contribs])]
            )

            self.log( "Highest absolute contribution of primitive internal coordinate "
                     f"{self.prim_coord} in mode {self.root:02d}.")
            self.log(f"Absolute contributions > {self.prim_contrib_thresh:.04f}"
                     f"\n{contrib_str}")
            used_str = f"contribution of primitive coordinate {self.prim_coord}"
        elif self.hessian_ref is not None:
            eigvals_ref, eigvecs_ref = np.linalg.eigh(self.hessian_ref)
            self.log_negative_eigenvalues(eigvals_ref, "Reference ")
            assert eigvals_ref[0] < -self.small_eigval_thresh
            ref_mode = eigvecs_ref[:,0]
            overlaps = np.einsum("ij,j->i", eigvecs.T, ref_mode)
            ovlp_str = np.array2string(overlaps, precision=4)
            self.log( "Overlaps between eigenvectors of current hessian "
                     f"TS mode from reference hessian:")
            self.log(f"\t{ovlp_str}")

            self.root = np.abs(overlaps).argmax()
            print( "Highest overlap between reference TS mode and "
                  f"eigenvector {self.root:02d}.")
            used_str = "overlap with reference TS mode"
        else:
            used_str = f"root={self.root}"
        self.log(f"Used {used_str} to select inital TS mode.")

        # Check if the selected mode (root) is a sensible choice.
        #
        # small_eigval_thresh is positive and we dont take the absolute value
        # of the eigenvalues. So we multiply small_eigval_thresh to get a
        # negative number.
        assert eigvals[self.root] < -self.small_eigval_thresh, \
             "Expected negative eigenvalue! Eigenvalue of selected TS-mode " \
            f"{self.root:02d} is above the the threshold of " \
            f"{self.small_eigval_thresh:.6e}!"

        # Select an initial TS-mode by root index. self.root may have been
        # modified by using a reference hessian.
        self.ts_mode = eigvecs[:,self.root]
        self.log(f"Using mode {self.root:02d} with eigenvalue "
                 f"{eigvals[self.root]:.6f} as TS mode.")
        self.log("")

    def update_ts_mode(self, eigvals, eigvecs):
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        assert neg_num >= 1, \
            "Need at least 1 negative eigenvalue for TS optimization."
        self.log_negative_eigenvalues(eigvals)

        # Select TS mode with biggest overlap to the previous TS mode
        self.log("Overlaps of previous TS mode with current imaginary mode(s):")
        ovlps = [np.abs(imag_mode.dot(self.ts_mode)) for imag_mode in eigvecs.T[:neg_num]]
        for i, ovlp in enumerate(ovlps):
            self.log(f"\t{i:02d}: {ovlp:.6f}")
        max_ovlp_ind = np.argmax(ovlps)
        max_ovlp = ovlps[max_ovlp_ind]
        self.log(f"Highest overlap: {max_ovlp:.6f}, mode {max_ovlp_ind}")
        self.log(f"Continuing with mode {max_ovlp_ind} as TS mode.")
        self.root = max_ovlp_ind
        self.ts_mode = eigvecs.T[max_ovlp_ind]
