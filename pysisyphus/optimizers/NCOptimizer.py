# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1498468
#     Bour, 2002

import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.helpers_pure import eigval_to_wavenumber 


class NCOptimizer(HessianOptimizer):

    def __init__(self, geometry, *args, freeze_modes=None, **kwargs):
        super().__init__(geometry, **kwargs)

        assert not self.is_cos and self.geometry.coord_type == "cart", \
            "NCOptimizer can't be used with ChainOfStates-methods and " \
            "coordinate systems beside cartesians ('coord_type: cart')."

        self.freeze_modes = freeze_modes
        
        self.sqrt_m = np.sqrt(self.geometry.masses_rep)

    def optimize(self):
        grad = self.geometry.gradient
        self.forces.append(-grad.copy())
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

        # eigvals_org, eigvecs_org = np.linalg.eigh(self.H)
        # grad_trans = eigvals_org.T.dot(grad)

        # Mass-weighted hessian
        H_mw = self.geometry.mass_weigh_hessian(self.H)
        # Project out translation and rotation
        if self.geometry.coords.size > 3:
            H_mw = self.geometry.eckart_projection(H_mw)
        eigvals, eigvecs = np.linalg.eigh(H_mw)
        # Drop translational/rotational modes as they will have
        # small (zero) eigenvalues. Zero the corresponding gradient entries.
        eigvals, eigvecs, small_inds = self.filter_small_eigvals(eigvals, eigvecs, mask=True)

        wavenumbers = eigval_to_wavenumber(eigvals)

        freeze_inds = np.zeros_like(wavenumbers, dtype=bool)
        if self.freeze_modes:
            freeze_inds = wavenumbers < self.freeze_modes
            eigvals = eigvals[~freeze_inds]
            eigvecs = eigvecs[:,~freeze_inds]
            self.log(f"{np.sum(freeze_inds)} normal modes will be frozen.")

        frozen_str = ["(frozen)" if frozen else "" for frozen in freeze_inds]
        wavenumber_str = "\n".join([
                            f"\t{i:> 3d}: {wn:> 8.2f} cm⁻¹ {frz}"
                            for i, (wn, frz) in enumerate(zip(wavenumbers, frozen_str))
        ])
        self.log("Frequencies:\n" + wavenumber_str)

        # Transform gradient to eigensystem of Hessian. We also have to use
        # the mass-weighted gradient.
        grad_q = eigvecs.T.dot(self.geometry.mw_gradient)

        mw_H_aug = self.get_augmented_hessian(eigvals, grad_q)
        # Discard eigenvector for now
        mw_step, eigval, nu, _ = self.solve_rfo(mw_H_aug, "min")
        # Transform back to original basis. Right now the step is still in
        # mass-weighted coordinates.
        mw_step = eigvecs @ mw_step
        # Un-massweigh step
        step = mw_step / self.sqrt_m

        step_norm = np.linalg.norm(step)

        if step_norm > self.trust_radius:
            step = step / step_norm * self.trust_radius

        if self.freeze_modes:
            # With frozen modes we only want to consider the gradient contributions
            # from non-frozen modes. So here we transform back the gradient and
            # un-weigh it.
            grad = eigvecs.dot(grad_q) * self.sqrt_m
            self.modified_forces.append(-grad)

        quadratic_prediction = step @ grad + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step
