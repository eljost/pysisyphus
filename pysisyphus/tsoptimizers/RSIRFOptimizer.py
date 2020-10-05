# [1] https://doi.org/10.1007/s002140050387
#     Bofill, 1998


import numpy as np

from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


class RSIRFOptimizer(TSHessianOptimizer):

    def optimize(self):
        energy, gradient, H, eigvals, eigvecs, resetted = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)

        # Transform gradient to eigensystem of hessian
        gradient_trans = eigvecs.T.dot(gradient)
        # Minimize energy along all modes, except the TS-mode
        min_indices = [i for i in range(gradient_trans.size) if i != self.root]
        # Get line search steps, if requested.
        ip_step_trans, ip_gradient_trans = self.step_and_grad_from_line_search(
            energy, gradient_trans, eigvecs, min_indices
        )
        ip_gradient = eigvecs.dot(ip_gradient_trans)
        ip_step = eigvecs.dot(ip_step_trans)

        self.log(
            "Using projection to construct image potential gradient "
            f"and hessian for root {self.root}."
        )
        trans_vec = eigvecs[:, self.root]
        # Projection matrix to construct g* and H*
        P = np.eye(self.geometry.coords.size) - 2 * np.outer(trans_vec, trans_vec)
        H_star = P.dot(H)
        eigvals_, eigvecs_ = np.linalg.eigh(H_star)
        # Neglect small eigenvalues
        eigvals_, eigvecs_ = self.filter_small_eigvals(eigvals_, eigvecs_)

        # grad_star = P.dot(gradient)
        grad_star = P.dot(ip_gradient)
        step = self.get_rs_step(eigvals_, eigvecs_, grad_star, name="RS-I-RFO")
        step += ip_step

        quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step
