from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer


class SteepestDescent(BacktrackingOptimizer):
    def __init__(self, geometry, alpha=0.1, **kwargs):
        super().__init__(geometry, alpha=alpha, **kwargs)

    def optimize(self):

        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.skip = self.backtrack(self.forces[-1], self.forces[-2])

        step = self.alpha * self.forces[-1]
        step = self.scale_by_max_step(step)
        return step
