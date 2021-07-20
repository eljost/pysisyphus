import numpy as np

from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.intcoords.PrimTypes import normalize_prim_input, normalize_prim_inputs
from pysisyphus.optimizers import poly_fit
from pysisyphus.optimizers.guess_hessians import ts_hessian
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.Optimizer import get_data_model, get_h5_group


class TSHessianOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    def __init__(
        self,
        geometry,
        root=0,
        hessian_ref=None,
        prim_coord=None,
        rx_coords=None,
        rx_mode=None,
        hessian_init="calc",
        hessian_update="bofill",
        hessian_recalc_reset=True,
        max_micro_cycles=50,
        trust_radius=0.3,
        augment_bonds=False,
        min_line_search=False,
        max_line_search=False,
        **kwargs,
    ):

        assert (
            hessian_update == "bofill"
        ), "Bofill update is recommended in a TS-optimization."

        super().__init__(
            geometry,
            hessian_init=hessian_init,
            hessian_update=hessian_update,
            trust_radius=trust_radius,
            hessian_recalc_reset=hessian_recalc_reset,
            **kwargs,
        )

        self.root = int(root)
        self.hessian_ref = hessian_ref
        try:
            self.hessian_ref = np.loadtxt(self.hessian_ref)
            _ = self.geometry.coords.size
            expected_shape = (_, _)
            shape = self.hessian_ref.shape
            assert shape == expected_shape, (
                f"Shape of reference Hessian {shape} doesn't match the expected "
                f"shape {expected_shape} of the Hessian for the current coordinates!"
            )
        except OSError as err:
            self.log(
                f"Tried to load reference Hessian from '{self.hessian_ref}' "
                "but the file could not be found."
            )
        except ValueError as err:
            self.log(f"No reference Hessian provided.")

        # Select initial root according to highest contribution of 'prim_coord'
        if prim_coord is not None:
            prim_coord = normalize_prim_input(prim_coord)[0]
        self.prim_coord = prim_coord

        # Construct initial imaginary mode according the primitive internals in
        # 'rx_coords' by modifying a model Hessian.
        if rx_coords is not None:
            rx_coords = normalize_prim_inputs(rx_coords)
        self.rx_coords = rx_coords

        # Construct initial imaginary mode from the given inputs while respecting
        # the given weight/phase factors.
        self.rx_mode = rx_mode

        # Bond augmentation is only useful with calculated hessians
        self.augment_bonds = augment_bonds and (self.hessian_init == "calc")
        self.min_line_search = min_line_search
        self.max_line_search = max_line_search

        self.ts_mode = None
        self.max_micro_cycles = max_micro_cycles
        self.prim_contrib_thresh = 0.05

        self.alpha0 = 1

    def prepare_opt(self, *args, **kwargs):
        if self.augment_bonds:
            self.geometry = augment_bonds(self.geometry, root=self.root)
            # Update data model and HD5 shapes, as the number of coordinates
            # may have changed.
            if self.dump:
                self.data_model = get_data_model(
                    self.geometry, self.is_cos, self.max_cycles
                )
                self.h5_group = get_h5_group(
                    self.h5_fn, self.h5_group_name, self.data_model
                )

        # Calculate/set initial hessian
        super().prepare_opt(*args, **kwargs)

        # Assume a guess hessian when not calculated. This hessian has to be
        # modified according to the assumed reaction coordinates.
        if self.hessian_init != "calc" and (self.rx_coords is not None):
            assert self.geometry.coord_type != "cart", (
                "Using a modified guess Hessian for TS-optimizations is "
                "only supported in redundand internal coordinates "
                "(coord_type=redund)"
            )
            prim_inds = [
                self.geometry.internal.get_index_of_typed_prim(rxc)
                for rxc in self.rx_coords
            ]
            missing_prim_inds = [
                self.rx_coords[i] for i, _ in enumerate(prim_inds) if _ is None
            ]
            assert len(missing_prim_inds) == 0, (
                "Some of the requested reaction coordinates are not defined: "
                f"{missing_prim_inds}"
            )
            self.H = ts_hessian(self.H, coord_inds=prim_inds)

        # Determiniation of initial mode either by using a provided
        # reference hessian, or by using a supplied root.

        eigvals, eigvecs = np.linalg.eigh(self.H)
        neg_inds = eigvals < -self.small_eigval_thresh
        self.log_negative_eigenvalues(eigvals, "Initial ")

        self.log("Determining initial TS mode to follow uphill.")
        # Select an initial TS-mode by highest overlap with eigenvectors from
        # reference Hessian.
        if self.prim_coord:
            prim_ind = self.geometry.internal.get_index_of_typed_prim(self.prim_coord)
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
                print(
                    "No imaginary mode with significant contribution "
                    f"(>{self.prim_contrib_thresh:.3f}) of primitive internal "
                    f"{self.prim_coord} found!"
                )
                raise err
            self.root = max_contrib_ind

            contrib_str = "\n".join(
                [
                    f"\t{ind:03d}: {contrib:.4f}"
                    for ind, contrib in zip(big_inds, np.abs(prim_row)[big_contribs])
                ]
            )

            self.log(
                "Highest absolute contribution of primitive internal coordinate "
                f"{self.prim_coord} in mode {self.root:02d}."
            )
            self.log(
                f"Absolute contributions > {self.prim_contrib_thresh:.04f}"
                f"\n{contrib_str}"
            )
            used_str = f"contribution of primitive coordinate {self.prim_coord}"
        elif self.hessian_ref is not None:
            eigvals_ref, eigvecs_ref = np.linalg.eigh(self.hessian_ref)
            self.log_negative_eigenvalues(eigvals_ref, "Reference ")
            assert eigvals_ref[0] < -self.small_eigval_thresh
            ref_mode = eigvecs_ref[:, 0]
            overlaps = np.einsum("ij,j->i", eigvecs.T, ref_mode)
            ovlp_str = np.array2string(overlaps, precision=4)
            self.log(
                "Overlaps between eigenvectors of current Hessian "
                f"TS mode from reference Hessian:"
            )
            self.log(f"\t{ovlp_str}")

            self.root = np.abs(overlaps).argmax()
            print(
                "Highest overlap between reference TS mode and "
                f"eigenvector {self.root:02d}."
            )
            used_str = "overlap with reference TS mode"
        # Construct an approximate initial mode according to user input
        # and calculate overlaps with the current eigenvectors.
        elif self.rx_mode is not None:
            self.log(f"Constructing reference mode, according to user input")
            assert self.geometry.coord_type != "cart"
            mode = np.zeros_like(self.geometry.coords)
            int_ = self.geometry.internal
            for typed_prim, val in self.rx_mode:
                typed_prim = normalize_prim_input(typed_prim)[0]
                ind = int_.get_index_of_typed_prim(typed_prim)
                mode[ind] = val
                self.log(f"\tIndex {ind} (coordinate {typed_prim}) set to {val:.4f}")
            mode /= np.linalg.norm(mode)
            mode_str = np.array2string(mode, precision=2)
            self.log(f"Final, normalized, reference mode: {mode_str}")

            # Calculate overlaps in non-redundant subspace by zeroing overlaps
            # in the redundant subspace.
            small_inds = np.abs(eigvals) < self.small_eigval_thresh
            # Take absolute value, because sign of eigenvectors is ambiguous.
            ovlps = np.abs(np.einsum("ij,i->j", eigvecs, mode))
            ovlps[small_inds] = 0.0
            self.root = np.argmax(ovlps)
            used_str = "overlap with user-generated mode"
        else:
            used_str = f"root={self.root}"
        self.log(f"Used {used_str} to select inital TS mode.")

        # This is currently commented out, as we can also start from
        # the convex region of the PES and maximize along a mode with
        # an initially positive eigenvalue.
        #
        # Check if the selected mode (root) is a sensible choice.
        #
        # small_eigval_thresh is positive and we dont take the absolute value
        # of the eigenvalues. So we multiply small_eigval_thresh to get a
        # negative number.
        # assert eigvals[self.root] < -self.small_eigval_thresh, (
        # "Expected negative eigenvalue! Eigenvalue of selected TS-mode "
        # f"{self.root:02d} is above the the threshold of "
        # f"{self.small_eigval_thresh:.6e}!"
        # )

        # Select an initial TS-mode by root index. self.root may have been
        # modified by using a reference hessian.
        self.ts_mode = eigvecs[:, self.root]
        self.ts_mode_eigval = eigvals[self.root]
        self.log(
            f"Using root {self.root:02d} with eigenvalue "
            f"{self.ts_mode_eigval:.6f} as TS mode.\n"
        )

    def update_ts_mode(self, eigvals, eigvecs):
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        self.log_negative_eigenvalues(eigvals)

        # When we left the convex region of the PES we only compare to other
        # imaginary modes ... is this a bad idea? Maybe we should use all modes
        # for the overlaps?!
        if self.ts_mode_eigval < 0:
            infix = "imaginary "
            ovlp_eigvecs = eigvecs[:, :neg_num]
            eigvals = eigvals[:neg_num]
            # When the eigenvalue corresponding to the TS mode has been negative once,
            # we should not lose all negative eigenvalues. If this happens something went
            # wrong and we crash :)
            assert (
                neg_num >= 1
            ), "Need at least 1 negative eigenvalue for TS optimization."
        # Use all eigenvectors for overlaps when the eigenvalue corresponding to the TS
        # mode is still positive.
        else:
            infix = ""
            ovlp_eigvecs = eigvecs

        # Select new TS mode according to biggest overlap with previous TS mode.
        self.log(f"Overlaps of previous TS mode with current {infix}mode(s):")
        ovlps = np.abs(np.einsum("ij,i->j", ovlp_eigvecs, self.ts_mode))
        for i, ovlp in enumerate(ovlps):
            self.log(f"\t{i:02d}: {ovlp:.6f}")
        max_ovlp_ind = np.argmax(ovlps)
        max_ovlp = ovlps[max_ovlp_ind]
        self.log(f"Highest overlap: {max_ovlp:.6f}, mode {max_ovlp_ind}")
        self.log(f"Maximizing along TS mode {max_ovlp_ind}.")
        self.root = max_ovlp_ind
        self.ts_mode = ovlp_eigvecs.T[self.root]
        self.ts_mode_eigval = eigvals[self.root]

    @staticmethod
    def do_line_search(e0, e1, g0, g1, prev_step, maximize, logger=None):
        poly_fit_kwargs = {
            "e0": e0,
            "e1": e1,
            "g0": g0,
            "g1": g1,
            "maximize": maximize,
        }
        if not maximize:
            poly_fit_kwargs.update(
                {
                    "g0": prev_step.dot(g0),
                    "g1": prev_step.dot(g1),
                }
            )
        prefix = "Max" if maximize else "Min"

        fit_result = poly_fit.quartic_fit(**poly_fit_kwargs)
        fit_energy = None
        fit_grad = None
        fit_step = None
        if fit_result and (0.0 < fit_result.x <= 2.0):
            x = fit_result.x
            log(logger, f"{prefix}-subpsace interpolation succeeded: x={x:.6f}")
            fit_energy = fit_result.y
            fit_step = (1 - x) * -prev_step
            fit_grad = (1 - x) * g0 + x * g1
        return fit_energy, fit_grad, fit_step

    def step_and_grad_from_line_search(
        self, energy, gradient_trans, eigvecs, min_indices
    ):
        ip_step = np.zeros_like(gradient_trans)
        ip_gradient_trans = gradient_trans.copy()

        if self.max_line_search and self.cur_cycle > 0:
            prev_energy = self.energies[-2]
            prev_gradient = -self.forces[-2]
            prev_gradient_trans = eigvecs.T.dot(prev_gradient)
            prev_step = self.steps[-1]
            prev_step_trans = eigvecs.T.dot(prev_step)

            # Max subspace
            # max_energy, max_gradient, max_step = self.do_max_line_search(
            max_energy, max_gradient, max_step = self.do_line_search(
                prev_energy,
                energy,
                prev_gradient_trans[self.root],
                gradient_trans[self.root],
                prev_step=prev_step_trans[self.root],
                maximize=True,
                logger=self.logger,
            )
            if max_gradient is not None:
                ip_gradient_trans[self.root] = max_gradient
                ip_step[self.root] = max_step

        if self.min_line_search and self.cur_cycle > 0:
            # Min subspace
            # min_energy, min_gradient, min_step = self.do_min_line_search(
            min_energy, min_gradient, min_step = self.do_line_search(
                prev_energy,
                energy,
                prev_gradient_trans[min_indices],
                gradient_trans[min_indices],
                prev_step=prev_step_trans[min_indices],
                maximize=False,
                logger=self.logger,
            )
            if min_gradient is not None:
                ip_gradient_trans[min_indices] = min_gradient
                ip_step[min_indices] = min_step
        return ip_step, ip_gradient_trans
