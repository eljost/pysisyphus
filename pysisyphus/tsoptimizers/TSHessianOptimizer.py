from typing import List, Optional

import h5py
import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import log
from pysisyphus.intcoords.augment_bonds import augment_bonds
from pysisyphus.intcoords.PrimTypes import normalize_prim_input, normalize_prim_inputs
from pysisyphus.optimizers import poly_fit
from pysisyphus.optimizers.guess_hessians import ts_hessian, HessInit
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer, HessUpdate
from pysisyphus.optimizers.Optimizer import get_data_model, get_h5_group


class TSHessianOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    valid_updates = ("bofill", "ts_bfgs", "ts_bfgs_org", "ts_bfgs_rev")

    def __init__(
        self,
        geometry: Geometry,
        roots: Optional[List[int]] = None,
        root: int = 0,
        hessian_ref: Optional[str] = None,
        rx_modes=None,
        prim_coord=None,
        rx_coords=None,
        hessian_init: HessInit = "calc",
        hessian_update: HessUpdate = "bofill",
        hessian_recalc_reset: bool = True,
        max_micro_cycles: int = 50,
        trust_radius: float = 0.3,
        trust_max: float = 0.5,
        augment_bonds: bool = False,
        min_line_search: bool = False,
        max_line_search: bool = False,
        assert_neg_eigval: bool = False,
        **kwargs,
    ) -> None:
        """Baseclass for transition state optimizers utilizing Hessian information.

        Several arguments expect a typed primitive or an iterable of typed primitives.
        A typed primitive is specified as (PrimType, int, int, ...), e.g., for a bond
        between atoms 0 and 1: (BOND, 0, 1) or for a bend between the atom triple 0, 1, 2
        as (BEND, 0, 1, 2).

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        roots
            Indices of modes to maximize along, e.g., to optimize saddle points of 2nd order.
            Overrides 'root'.
        root
            Index of imaginary mode to maximize along. Shortcut for 'roots' with only one root.
        hessian_ref
            Filename pointing to a pysisyphus HDF5 Hessian.
        rx_modes : iterable of (iterable of (typed_prim, phase_factor))
            Select initial root(s) by overlap with a modes constructed from the given
            typed primitives with respective phase factors.
        prim_coord : typed_prim
            Select initial root/imaginary mode by overlap with this internal coordinate.
            Shortcut for 'rx_modes' with only one internal coordinate.
        rx_coords : iterable of (typed_prim)
            Construct imaginary mode comprising the given typed prims by modifying
            a model Hessian.
        hessian_init
            Type of initial model Hessian.
        hessian_update
            Type of Hessian update. Defaults to BFGS for minimizations and Bofill
            for saddle point searches.
        hessian_recalc_reset
            Whether to skip Hessian recalculation after reset. Undocumented.
        max_micro_cycles
            Maximum number of RS iterations.
        trust_radius
            Initial trust radius in whatever unit the optimization is carried out.
        trust_max
            Maximum trust radius.
        augment_bonds
            Try to derive additional streching coordinates from the imaginary mode.
        min_line_search
            Carry out line search along the imaginary mode.
        max_line_search
            Carry out line search in the subspace that is minimized.
        assert_neg_eigval
            Check for the existences for at least one significant negative eigenvalue.
            If enabled and no negative eigenvalue is present the optimization will be
            aborted.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the HessianOptimizer/Optimizer baseclass.
        """

        assert (
            hessian_update in self.valid_updates
        ), f"Invalid Hessian update. Please chose from: {self.valid_updates}!"

        super().__init__(
            geometry,
            hessian_init=hessian_init,
            hessian_update=hessian_update,
            trust_radius=trust_radius,
            trust_max=trust_max,
            hessian_recalc_reset=hessian_recalc_reset,
            **kwargs,
        )

        if (root is not None) and (roots is None):
            roots = [
                root,
            ]
        elif roots is None:
            roots = list()
        self.roots = roots
        self.log(f"{self.roots=}")
        self.hessian_ref = hessian_ref
        try:
            with h5py.File(self.hessian_ref, "r") as handle:
                self.hessian_ref = handle["hessian"][:]
            _ = self.geometry.coords.size
            expected_shape = (_, _)
            shape = self.hessian_ref.shape
            # Hessian is not yet converted to the correct coordinate system if
            # coord_type != cart.
            assert (
                self.geometry.coord_type == "cart"
            ), "hessian_ref with internal coordinates are not yet supported."
            assert shape == expected_shape, (
                f"Shape of reference Hessian {shape} doesn't match the expected "
                f"shape {expected_shape} of the Hessian for the current coordinates!"
            )
        except OSError:
            self.log(
                f"Tried to load reference Hessian from '{self.hessian_ref}' "
                "but the file could not be found."
            )
            self.hessian_ref = None
        except (ValueError, TypeError):
            self.log("No reference Hessian provided.")

        # Select initial root according to highest contribution of 'prim_coord'
        if prim_coord is not None:
            self.log("'prim_coord' given. Overwriting/setting 'rx_modes'.")
            rx_modes = [[[prim_coord, 1.0]]]
        self.prim_coord = prim_coord

        # Construct initial imaginary mode from the given inputs while respecting
        # the given weight/phase factors.
        self.rx_modes = rx_modes

        # Construct initial imaginary mode according the primitive internals in
        # 'rx_coords' by modifying a model Hessian.
        if rx_coords is not None:
            rx_coords = normalize_prim_inputs(rx_coords)
        self.rx_coords = rx_coords

        # Bond augmentation is only useful with calculated hessians
        self.augment_bonds = augment_bonds and (self.hessian_init == "calc")
        self.min_line_search = min_line_search
        self.max_line_search = max_line_search
        self.assert_neg_eigval = assert_neg_eigval

        self.ts_modes = list()
        self.max_micro_cycles = max_micro_cycles
        self.prim_contrib_thresh = 0.05

        self.alpha0 = 1

    @property
    def root(self):
        return self.roots[0]

    @root.setter
    def root(self, root):
        raise Exception("Setting 'self.root' is deprecated. Set 'self.roots' instead.")

    @property
    def roots(self):
        return self._roots

    @roots.setter
    def roots(self, roots):
        self._roots = np.array(roots, dtype=int)

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
        if self.hessian_ref is not None:
            eigvals_ref, eigvecs_ref = np.linalg.eigh(self.hessian_ref)
            self.log_negative_eigenvalues(eigvals_ref, "Reference ")
            assert eigvals_ref[0] < -self.small_eigval_thresh
            ref_mode = eigvecs_ref[:, 0]
            overlaps = np.einsum("ij,j->i", eigvecs.T, ref_mode)
            ovlp_str = np.array2string(overlaps, precision=4)
            self.log(
                "Overlaps between eigenvectors of current Hessian "
                "TS mode from reference Hessian:"
            )
            self.log(f"\t{ovlp_str}")

            root = np.abs(overlaps).argmax()
            self.roots = [
                root,
            ]
            print(
                "Highest overlap between reference TS mode and "
                f"eigenvector {self.root:02d}."
            )
            used_str = "overlap with reference TS mode"
        # Construct an approximate initial mode according to user input
        # and calculate overlaps with the current eigenvectors.
        elif self.rx_modes is not None:
            self.log("Constructing reference mode, according to user input")
            assert self.geometry.coord_type != "cart"
            int_ = self.geometry.internal
            modes = list()
            for i, rx_mode in enumerate(self.rx_modes):
                mode = np.zeros_like(self.geometry.coords)
                for typed_prim, val in rx_mode:
                    typed_prim = normalize_prim_input(typed_prim)[0]
                    ind = int_.get_index_of_typed_prim(typed_prim)
                    mode[ind] = val
                    self.log(
                        f"\tIndex {ind} (coordinate {typed_prim}) set to {val:.4f}"
                    )
                mode /= np.linalg.norm(mode)
                modes.append(mode)
                mode_str = np.array2string(mode, precision=2)
                self.log(f"Normalized reference mode {i:02d}: {mode_str}")

            # Calculate overlaps in non-redundant subspace by zeroing overlaps
            # in the redundant subspace.
            small_inds = np.abs(eigvals) < self.small_eigval_thresh
            # Take absolute value, because sign of eigenvectors is ambiguous.
            ovlps = np.abs(np.einsum("ij,ki->kj", eigvecs, modes))
            ovlps[:, small_inds] = 0.0
            self.roots = ovlps.argmax(axis=1)
            used_str = "overlap with user-generated mode"
        else:
            used_str = f"root(s)={self.roots}"
        self.log(f"Used {used_str} to select inital TS mode.")

        # Below, some code is found, that checks if the chosen root(s) are a
        # sensible choice, i.e., if they are negative. Currently, it is commented out,
        # as we can also start from the convex region of the PES.
        #
        # Check if the selected mode (root) is a sensible choice.
        #
        # small_eigval_thresh is positive and we dont take the absolute value
        # of the eigenvalues. So we multiply small_eigval_thresh to get a
        # negative number.
        # assert (eigvals[self.roots] < -self.small_eigval_thresh).all(), (
        # "Expected negative eigenvalue(s)! Eigenvalues of selected TS-modes "
        # f"are above the the threshold of self.small_eigval_thresh:.6e}!"
        # )

        # Select an initial TS-mode by root index. self.root may have been
        # modified by using a reference hessian.
        self.ts_modes = eigvecs[:, self.roots].T
        self.ts_mode_eigvals = eigvals[self.roots]
        self.log(
            f"Using root(s) {self.roots} with eigenvalues "
            f"{np.array2string(self.ts_mode_eigvals, precision=6)} as TS mode.\n"
        )

    def update_ts_mode(self, eigvals, eigvecs):
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        self.log_negative_eigenvalues(eigvals)

        # When we left the convex region of the PES we only compare to other
        # imaginary modes ... is this a bad idea? Maybe we should use all modes
        # for the overlaps?!
        if (self.ts_mode_eigvals < 0).all() and neg_num > 0:
            infix = "imaginary "
            ovlp_eigvecs = eigvecs[:, :neg_num]
            eigvals = eigvals[:neg_num]
        # When the eigenvalue corresponding to the TS mode has been negative once,
        # we should not lose all negative eigenvalues. If this happens something went
        # wrong and we crash :)
        elif self.assert_neg_eigval and neg_num == 0:
            raise AssertionError(
                "Need at least 1 negative eigenvalue for TS optimization."
            )
        # Use all eigenvectors for overlaps when the eigenvalue corresponding to the TS
        # mode is still positive.
        else:
            infix = ""
            ovlp_eigvecs = eigvecs

        # Select new TS mode according to biggest overlap with previous TS mode.
        self.log(f"Overlaps of previous TS mode with current {infix}mode(s):")
        ovlps = np.abs(np.einsum("ij,ki->kj", ovlp_eigvecs, self.ts_modes))
        for i, ovlp in enumerate(ovlps):
            self.log(f"\tTS mode {i:02d}: {np.array2string(ovlp, precision=3)}")
        max_ovlp_inds = np.argmax(ovlps, axis=1)
        for i, _ in enumerate(self.ts_modes):
            max_ovlp_ind = max_ovlp_inds[i].argmax()
            self.log(
                f"Mode {i}: highest overlap: {ovlps[i, max_ovlp_ind]:.6f} with mode "
                f"{max_ovlp_ind}"
            )
        self.roots = max_ovlp_inds
        self.ts_modes = ovlp_eigvecs.T[self.roots]
        self.ts_mode_eigvals = eigvals[self.roots]

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
        self,
        energy,
        gradient_trans,
        eigvecs,
        min_indices,
        max_indices,
    ):
        ip_step = np.zeros_like(gradient_trans)
        ip_gradient_trans = gradient_trans.copy()

        if self.max_line_search and self.cur_cycle > 0:
            prev_energy = self.energies[-2]
            prev_gradient = -self.forces[-2]
            try:
                prev_gradient_trans = eigvecs.T.dot(prev_gradient)
                prev_step = self.steps[-1]
                prev_step_trans = eigvecs.T.dot(prev_step)
            # Will be raised when coordinates were rebuilt and the array shapes differe.
            except ValueError:
                return ip_step, ip_gradient_trans

            # Max subspace
            # max_energy, max_gradient, max_step = self.do_max_line_search(
            max_energy, max_gradient, max_step = self.do_line_search(
                prev_energy,
                energy,
                prev_gradient_trans[max_indices],
                gradient_trans[max_indices],
                prev_step=prev_step_trans[max_indices],
                maximize=True,
                logger=self.logger,
            )
            if max_gradient is not None:
                ip_gradient_trans[max_indices] = max_gradient
                ip_step[max_indices] = max_step

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
