import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class GAD(HessianOptimizer):

    def __init__(self, geom, dt=0.1, **kwargs):
        super().__init__(geom, **kwargs)

        self.dt = dt

    def prepare_opt(self):
        super().prepare_opt()

        # Usually eigenvector of H matrix
        v = self.geometry.gradient
        self.v = v / np.linalg.norm(v)

        self.I = np.eye(self.geometry.coords.size)

    def optimize(self):
        grad = self.geometry.gradient
        self.forces.append(-grad)
        hessian = self.geometry.hessian

        P = np.outer(self.v, self.v)

        dxdt = -(self.I - 2*P) @ grad
        dvdt = -(self.I - P) @ hessian @ self.v
        dv = dvdt * self.dt
        v = self.v + dv
        self.v = v / np.linalg.norm(v)

        step = dxdt * self.dt
        step_norm = np.linalg.norm(step)
        if step_norm > self.dt:
            step = self.dt * step / step_norm
        return step
