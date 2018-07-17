#!/usr/bin/env python3

from pysisyphus.cos.NEB import NEB


class FreeEndNEB(NEB):
    pass

    @NEB.forces.getter
    def forces(self):
        forces = super().forces
        forces_last = self.images[-1].forces
        par_forces_last = self.get_parallel_forces(len(self.images)-1)
        lambda_ = (forces_last.dot(forces_last)
                   / -forces_last.dot(par_forces_last))
        fe_forces = forces_last + lambda_ * par_forces_last

        lambda_2 = (par_forces_last.dot(-forces_last) / forces_last.dot(forces_last))
        fe_forces_2 = par_forces_last + lambda_2*forces_last
        # self._forces = np.array(total_forces).flatten()
        self._forces[-forces_last.shape[0]:] = fe_forces_2
        self._forces[-forces_last.shape[0]:] = fe_forces

        # import pdb; pdb.set_trace()
        # print("huhu")
        return forces
