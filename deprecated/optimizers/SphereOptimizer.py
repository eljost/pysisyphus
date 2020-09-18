#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize

from pysisyphus.optimizers.Optimizer import Optimizer


class SphereOptimizer(Optimizer):

    def __init__(self, geometry, radius=0.3, **kwargs):
        super().__init__(geometry, **kwargs)

        self.radius = radius
        # self.initial_coords = self.geometry.coords.copy()
        self.initial_coords = np.array((0.05, 0, 0))
        self.ref_ind = 1

    def optimize(self):
        forces = self.geometry.forces
        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        ref_force = forces[self.ref_ind]
        cur_coords = self.geometry.coords
        ref_diff = cur_coords[self.ref_ind] - self.initial_coords[self.ref_ind]
        constr_forces = (
            forces - ref_force*(cur_coords - self.initial_coords) / ref_diff
        )
        direction = constr_forces / np.linalg.norm(constr_forces)

        step = 0.01*direction

        return step


def run():
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from pysisyphus.calculators.CerjanMiller import CerjanMiller
    from pysisyphus.calculators.AnaPot import AnaPot

    ref_coords = np.array((0.05, 0, 0))
    geom = CerjanMiller.get_geom(ref_coords)
    # geom = AnaPot.get_geom((-1, 1, 0))

    ref_ind = 1
    free_inds = np.arange(geom.coords.size) != ref_ind
    free_ref_coords = ref_coords[free_inds]
    ref_coord = ref_coords[ref_ind]

    def get_constr_coord(coords, radius):
        free_coords = coords[free_inds]
        _ = max(0, radius**2 - radius*np.sum((free_coords - free_ref_coords)**2))
        sqrt = np.sqrt(_)
        return ref_coord + sqrt

    def get_forces_mod(forces, coords):
        ref_force = forces[ref_ind]
        constr_coord = coords[ref_ind]
        return forces - ref_force*(coords-ref_coords)/(constr_coord-ref_coord)

    dR = 0.01
    alpha = 0.01
    R = dR
    # def wrapper(coords):
        # geom.coords = coords
        # constr_coord = get_constr_coord(geom.coords, R)
        # diff = R - np.linalg.norm(geom.coords - ref_coords)
        # geom.set_coord(ref_ind, constr_coord)
        # forces = geom.forces
        # forces_mod = get_forces_mod(forces, geom.coords)
        # rms = np.sqrt(np.mean(forces_mod**2))
        # return rms, -forces_mod

    # res = minimize(wrapper, geom.coords, jac=True)
    # print(res)
    # return

    dR = 0.01
    alpha = 0.01
    R = dR
    all_coords = list()
    radii = [R, ]
    for i in range(10):
        all_coords.append(geom.coords.copy())
        print(f"Cycle {i:02d}, R={R:.04f}")
        constr_coord = get_constr_coord(geom.coords, R)

        # Set constrained coordinate
        geom.set_coord(ref_ind, constr_coord)
        assert geom.coords[ref_ind] == constr_coord
        all_coords.append(geom.coords.copy())

        diff = R - np.linalg.norm(geom.coords - ref_coords)

        # Optimize remaining coordinates
        forces = geom.forces
        norm_forces = np.linalg.norm(forces)

        # Force acting on the constrained coordinate
        ref_force = forces[ref_ind]
        forces_mod_ = forces - ref_force*(geom.coords-ref_coords)/(constr_coord-ref_coord)
        forces_mod = get_forces_mod(forces, geom.coords)
        np.testing.assert_allclose(forces_mod_, forces_mod)
        norm_forces_mod = np.linalg.norm(forces_mod)
        import pdb; pdb.set_trace()

        print(f"\tnorm(f)={norm_forces:.6f} "
              f"norm(fm)={norm_forces_mod:.6f} dR={diff:.6f}"
        )
        # if norm_forces_mod < alpha:
            # # Increase radius
            # R += dR
            # radii.append(R)

            # # Start extrapolation if enough points are present
            # if len(radii) < 2:
                # continue
            # V = (all_coords[-1] - all_coords[-2]) / (radii[-1] - radii[-2])
            # import pdb; pdb.set_trace()
            # while True:
                # print("extrapol")
                # extrapolated_coords = geom.coords + V*dR
                # geom.coords = extrapolated_coords
                # forces = geom.forces
                # forces_norm = np.linalg.norm(forces)
                # forces_mod = get_forces_mod(forces, geom.coords)
                # norm_force_mod = np.linalg.norm(forces_mod)
                # if norm_forces_mod > alpha:
                    # break

        import pdb; pdb.set_trace()
        dir_ = forces_mod / np.linalg.norm(forces_mod)
        step = 0.01 * dir_
        new_coords = geom.coords + step
        geom.coords = new_coords

    pot = geom.calculator
    pot.plot()
    all_coords = np.array(all_coords)
    radii = np.array(radii)
    pot.ax.plot(all_coords[:,0], all_coords[:,1], "o-")
    pot.ax.scatter(ref_coords[0], ref_coords[1], s=20, c="r")
    for radius in radii:
        circle = Circle(ref_coords[:2], radius=radius, fill=False)
        pot.ax.add_artist(circle)
    for i, xy in enumerate(all_coords[:,:2]):
        pot.ax.annotate(i, xy)
    pot.ax.set_xlim(-0.2, 0.2)
    pot.ax.set_ylim(-0.2, 0.2)
    plt.show()


if __name__ == "__main__":
    run()
