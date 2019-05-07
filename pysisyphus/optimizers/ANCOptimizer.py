#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class ANCOptimizer(HessianOptimizer):

    def __init__(self, geometry, **kwargs):
        super().__init__(geometry, **kwargs)
        assert not self.is_cos and self.geometry.coord_type == "cart", \
            "ANCOptimizer can't be used with ChainOfStates-methods and " \
            "coordinate systems beside cartesians ('coord_type: cart')."

    def prepare_opt(self):
        super().prepare_opt()
        self.M_inv = self.geometry.mm_inv
    
    def optimize(self):
        grad = self.geometry.gradient
        self.forces.append(-grad.copy())
        self.energies.append(self.geometry.energy)

        H_mw = self.M_inv.dot(self.H).dot(self.M_inv)
        H_mw_proj = self.geometry.eckart_projection(H_mw)
        eigvals_mw, eigvecs_mw = np.linalg.eigh(H_mw_proj)

        # Neglect translational/rotational modes
        keep = np.abs(eigvals_mw) > 1e-12
        eigvals_mw = eigvals_mw[keep]
        eigvecs_mw = eigvecs_mw[:,keep]

        # Unweight mass-weighted normal modes
        eigvecs = self.M_inv.dot(eigvecs_mw)
        # Transform gradient
        grad_q = eigvecs.T.dot(grad)

        if self.cur_cycle > 0:
            self.update_hessian()

        dQ = -2*grad_q / (eigvals_mw + np.sqrt(eigvals_mw**2 + 4*(grad_q**2)))

        step = eigvecs.dot(dQ)
        step_norm = np.linalg.norm(step)

        if step_norm > self.trust_radius:
            step = step / step_norm * self.trust_radius
        return step
