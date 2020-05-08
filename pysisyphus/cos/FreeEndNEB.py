#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.NEB import NEB

# [1] https://www.pnas.org/content/pnas/104/9/3031.full.pdf
#     Zhu, 2006
#     Original method
# [2] http://dx.doi.org/10.1063/1.4962019
#     Zhang, 2016
#     FreeEnd Adaptive NEB


class FreeEndNEB(NEB):
    kinds = ("original", "improved", "improved_v2")

    def __init__(self, *args, kind="improved", fix_first=False, fix_last=False,
                 **kwargs):
        super().__init__(*args, fix_first=fix_first, fix_last=fix_last, **kwargs)

        assert self.fix_ends == False
        assert kind in self.kinds

        self.kind = kind

    @NEB.forces.getter
    def forces(self):
        forces = super().forces
        i_last = self.last_index
        forces_last = self.images[-1].forces
        forces_size = forces_last.size

        tangent = self.get_tangent(i_last)
        coord_diff = self.images[i_last].coords - self.images[i_last-1].coords
        distance = np.linalg.norm(coord_diff)
        spring_constant = self.k[i_last-1]
        par_forces = -spring_constant * distance * tangent

        # Original FE-NEB forces on last image, Eq. (1) in [1]
        if self.kind == "original":
            end_forces = (par_forces
                          - par_forces.dot(forces_last) / np.linalg.norm(forces_last)
                          * forces_last
            )
        # Improved FE-NEB forces proposed by Zhang in [2], Eq. (7)
        elif self.kind == "improved":
            end_forces = forces_last - forces_last.dot(tangent) * tangent
            # a = self.get_perpendicular_forces(i_last)
            # import pdb; pdb.set_trace()
            # np.testing.assert_allclose(a, end_forces)
        # Another proposed FE-NEB force by Zhang in [2], Appendix B
        elif self.kind == "improved_v2":
            lambda_ = np.linalg.norm(forces_last)/-forces_last.dot(par_forces)
            self.log(f"λ={lambda_:.6f}")
            end_forces = forces_last + lambda_ * par_forces
        else:
            raise Exception(f"Invalid kind={self.kind}! Valid kinds are: {self.kinds}.")
        forces[-forces_size:] = end_forces
        self._forces = forces

        return self._forces
