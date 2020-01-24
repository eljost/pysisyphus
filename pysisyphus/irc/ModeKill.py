import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.helpers import eigval_to_wavenumber, do_final_hessian
from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import bofill_update
from pysisyphus.TableFormatter import TableFormatter


class ModeKill(IRC):

    def __init__(self, geometry, kill_inds, nu_thresh=-5.,
                 do_hess=True, hessian_update=True, **kwargs):
        # Use tight convergence criteria so the IRC/ModeKill doesn't
        # ent prematurely.
        super().__init__(geometry, downhill=True, rms_grad_thresh=1e-6,
                         **kwargs)

        assert self.geometry.coord_type == "cart"

        self.kill_inds = np.array(kill_inds, dtype=int)
        self.nu_thresh = float(nu_thresh)
        self.do_hess = do_hess
        self.hessian_update = hessian_update

        self.fmt = {"float": lambda f: f"{f:.4f}",}
        self.ovlp_thresh = .3
        self.indices = np.arange(self.geometry.coords.size)
        self.neg_nus = list()

    def prepare(self, *args, **kwargs):
        super().prepare(*args, **kwargs)

        self.mw_H = self.mw_hessian
        self.update_mw_down_step()

    def update_mw_down_step(self):
        w, v = np.linalg.eigh(self.mw_H)

        self.kill_modes = v[:,self.kill_inds]
        nus = eigval_to_wavenumber(w)
        assert all(nus[self.kill_inds] < self.nu_thresh), \
              "ModeKill is intended for removal of imaginary frequencies " \
             f"below {self.nu_thresh} cm⁻¹! The specified indices " \
             f"{self.kill_inds} contain modes with positive frequencies " \
             f"({nus[self.kill_inds]} cm⁻¹). Please choose different kill_inds!"

        """After diagonalization of the mass-weighted hessian the signs
        of the eigenvectors are arbitrary.
        We determine the correct sign of the eigenvector(s) from its
        overlap(s) with the gradient.
        To decrease the overall energy we have to step against the
        gradient, so the overlap of gradient and the respective eigen-
        vector(s) must be negative.
        If an overlap is positive we flip the sign of the corresponding
        eigenvector.
        """
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

    def step(self):
        if self.hessian_update and (self.cur_cycle > 0):
            # Hessian update with mass-weighted values
            dx = self.irc_mw_coords[-1] - self.irc_mw_coords[-2]
            dg = (self.irc_mw_gradients[-1] - self.irc_mw_gradients[-2])
            d_mw_H, key = bofill_update(self.mw_H, dx, dg)
            self.mw_H += d_mw_H

            # norm(dx) is probably self.step_length ;)
            norm_dx = np.linalg.norm(dx)
            norm_dg = np.linalg.norm(dg)
            self.log(f"Did {key} hessian update: norm(dx)={norm_dx:.4e}, "
                     f"norm(dg)={norm_dg:.4e}."
            )
        else:
            # Recalculate exact hessian
            self.mw_H = self.mw_hessian
            self.log("Recalculated exact hessian.")

        w, v = np.linalg.eigh(self.mw_H)
        # Overlaps between current normal modes and the modes we want to
        # remove.
        overlaps = np.abs(np.einsum("ij,ik->jk", self.kill_modes, v))
        self.log(f"Overlaps between original modes and current modes:")
        # overlaps contains one row per mode to remove
        for i, ovlps in enumerate(overlaps):
            above_thresh = ovlps > self.ovlp_thresh
            ovlp_str = " ".join(
                [f"{i:02d}: {o:.4f}" for i, o in zip(self.indices[above_thresh],
                                                     ovlps[above_thresh])
                ]
            )
            self.log(f"Org. mode {i:02d}:\n\t\t\t{ovlp_str}")

        nus = eigval_to_wavenumber(w)
        neg_inds = nus <= self.nu_thresh
        neg_nus = nus[neg_inds]
        self.neg_nus.append(neg_nus)
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
            self.log(f"Eigenvalue(s) of mode(s) {flipped} became positive!")
        self.converged = len(self.kill_inds) == 0

        if not self.converged:
            self.mw_coords += self.mw_down_step

    def postprocess(self):
        super().postprocess()

        if self.do_hess:
            print()
            do_final_hessian(self.geometry)

    def get_additional_print(self):
        neg_nus = np.array2string(self.neg_nus[-1], precision=2)
        return f"\timag. ῦ: {neg_nus} cm⁻¹"
