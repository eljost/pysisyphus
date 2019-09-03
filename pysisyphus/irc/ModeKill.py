import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.helpers import eigval_to_wavenumber
from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter


class ModeKill(IRC):

    def __init__(self, geometry, kill_inds, nu_thresh=-5., **kwargs):
        super().__init__(geometry, downhill=True, **kwargs)

        assert self.geometry.coord_type == "cart"

        self.kill_inds = np.array(kill_inds, dtype=int)
        self.nu_thresh = float(nu_thresh)

        self.fmt = {"float": lambda f: f"{f:.4f}",}
        self.ovlp_thresh = .3
        self.indices = np.arange(self.geometry.coords.size)

    def update_mw_down_step(self):
        w, v = np.linalg.eigh(self.mw_hessian)

        self.kill_modes = v[:,self.kill_inds]
        nus = eigval_to_wavenumber(w)
        assert all(nus[self.kill_inds] < self.nu_thresh), \
              "ModeKill is intended for removal of imaginary frequencies " \
             f"below {self.nu_thresh} cm⁻¹! But the specified indices " \
             f"{self.kill_inds} contain modes with positive frequencies " \
             f"({nus[self.kill_inds]} cm⁻¹)."
        # We determine the correct sign of the eigenvector(s) from its
        # overlap(s) with the gradient.
        # To decrease the overall energy we have to step against the
        # gradient, so the overlap of gradient and the respective eigen-
        # vector(s) must be negative.
        # If an overlap is positive we flip the sign of the corresponding
        # eigenvector.
        mw_grad = self.mw_gradient
        mw_grad_normed = mw_grad / np.linalg.norm(mw_grad)
        overlaps = np.einsum("ij,i->j", self.kill_modes, mw_grad_normed)
        self.log("Overlaps between gradient and eigenvectors:")
        self.log(overlaps)
        flip = overlaps > 0
        self.log("Eigenvector signs to be flipped:")
        self.log(str(flip))
        self.kill_modes[:,flip] *= -1
        # Create the step as the sum of the downhill steps along the modes
        # to remove.
        self.mw_down_step = (self.step_length*self.kill_modes).sum(axis=1)

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        self.update_mw_down_step()

    def step(self):
        grad = self.mw_gradient
        energy = self.energy
        mw_hessian = self.mw_hessian

        w, v = np.linalg.eigh(mw_hessian)
        # Overlaps between current normal modes and the modes we want to
        # remove.
        overlaps = np.abs(np.einsum("ij,ik->jk", self.kill_modes, v))
        self.log(f"Overlaps between original modes and current modes:")
        # overlaps contains one row per mode to remove
        for i, ovlps in enumerate(overlaps):
            above_thresh = ovlps > self.ovlp_thresh
            # ovlp_str = np.array2string(ovlp[:10], prefix=" "*18)
            ovlp_str = " ".join(
                [f"{i:02d}: {o:.4f}" for i, o in zip(self.indices[above_thresh],
                                                     ovlps[above_thresh])
                ]
            )
            self.log(f"\tOrg. mode {i:02d}: {ovlp_str}")

        nus = eigval_to_wavenumber(w)
        neg_inds = nus <= self.nu_thresh
        neg_nus = nus[neg_inds]
        self.log(f"Wavenumbers of imaginary modes (<= {self.nu_thresh} cm⁻¹):")
        self.log(f"{neg_nus} cm⁻¹")

        # Check if any eigenvalues became positive. If so remove them and update
        # the step. If no mode to remove is left we are finished and can signal
        # convergence.
        argmax = overlaps.argmax(axis=1)
        eigvals = w[argmax]
        pos_eigvals = eigvals > 0
        if any(pos_eigvals):
            # Only keep negative eigenvalues.
            flipped = self.kill_inds[pos_eigvals]
            self.kill_inds = self.kill_inds[~pos_eigvals]
            self.update_mw_down_step()
            self.log("Eigenvalue(s) of mode(s) {flipped} became positive!")
        self.converged = len(self.kill_inds) == 0

        if not self.converged:
            self.mw_coords += self.mw_down_step
