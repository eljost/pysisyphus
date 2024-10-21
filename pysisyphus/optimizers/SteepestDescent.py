from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer


class SteepestDescent(BacktrackingOptimizer):
    def __init__(self, geometry, alpha=0.1, **kwargs):
        super().__init__(geometry, alpha=alpha, **kwargs)

    def get_step(
        self, energy, forces, hessian=None, eigvals=None, eigvecs=None, resetted=None
    ):
        if self.cur_cycle > 0:
            self.skip = self.backtrack(forces, self.forces[-2])

        step = self.alpha * forces
        step = self.scale_by_max_step(step)
        return step
