import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.helpers import eigval_to_wavenumber
from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter


class ModeKill(IRC):

    def __init__(self, geometry, kill_inds=None, nu_thresh=-5., **kwargs):
        super().__init__(geometry, downhill=True, **kwargs)

        self.kill_inds = list(kill_inds) if kill_inds else list()
        self.nu_thresh = float(nu_thresh)

    def update_mw_down_step(self):
        w, v = np.linalg.eigh(self.mw_hessian)

        self.kill_modes = v[:,self.kill_inds]
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
        flip = overlaps > 0
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
        self.log(f"Overlaps between current normal modes and modes to remove:")
        self.log(str(overlaps[:,:6]))

        nus = eigval_to_wavenumber(w)
        neg_inds = nus < self.nu_thresh
        neg_nus = nus[neg_inds]
        self.log("Wavenumber of imaginary modes:")
        self.log(str(neg_nus))

        self.mw_coords += self.mw_down_step
