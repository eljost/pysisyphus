#!/usr/bin/env python3

from pysisyphus.cos.NEB import NEB

# [1] http://dx.doi.org/10.1063/1.4962019


class FreeEndNEB(NEB):
    def __init__(self, *args, mod=False, **kwargs):
        super().__init__(*args, **kwargs)

        self.mod = mod

    def mod_end_forces(self, i):
        """Equation (7) in [1]."""
        assert (i == 0) or (i == self.last_index)
        true_forces = self.images[i].forces
        tangent = self.get_tangent(i)
        return true_forces - true_forces.dot(tangent) * tangent

    @NEB.forces.getter
    def forces(self):
        forces = super().forces
        i_last = self.last_index
        forces_last = self.images[-1].forces
        forces_size = forces_last.size
        if self.mod:
            start_forces = self.mod_end_forces(0)
            end_forces = self.mod_end_forces(self.last_index)
            forces[:forces_size] = start_forces
            forces[-forces_size:] = end_forces
        else:
            par_forces_last = self.get_parallel_forces(i_last)
            lambda_ = -forces_last.dot(forces_last)/par_forces_last.dot(forces_last)
            fe_forces = lambda_ * par_forces_last + forces_last
            forces[-forces_size:] = fe_forces
        self._forces = forces

        return self._forces
