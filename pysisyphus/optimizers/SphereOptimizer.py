#!/usr/bin/env python3

import numpy as np

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
    geom = CerjanMiller.get_geom((0.05, 0, 0))
    import pdb; pdb.set_trace()

    ref_ind = 1
    ref_coords = np.array((0.05, 0, 0 ))
    free_inds = np.arange(geom.coords.size) != ref_ind
    free_ref_coords = ref_coords[free_inds]
    ref_coord = ref_coords[ref_ind]

    def get_constr_coord(coords, radius):
        free_coords = coords[free_inds]
        _ = np.sqrt(radius**2 - radius*np.sum((free_coords - free_ref_coords)**2))
        return float(ref_coord + _)

    # for radius in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
    c = list()
    rads = np.array((0.1, 0.2))
    rs = list()
    dR = 0.02
    alpha = 0.01
    R = dR
    for i in range(15):
        c.append(geom.coords.copy())
        print(f"Cycle {i:02d}, R={R:.04f}")
        constr_coord = get_constr_coord(geom.coords, R)
        # Set constrained coordinate
        new_coords = geom.coords.copy()
        new_coords[ref_ind] = constr_coord
        geom.coords = new_coords
        print("\t", geom.coords)
        geom.clear()

        # Optimize remaining coordinates
        forces = geom.forces
        norm_forces = np.linalg.norm(forces)
        diff = R - np.linalg.norm(geom.coords - ref_coords)
        free_forces = forces[free_inds]
        norm_free_forces = np.linalg.norm(free_forces)
        # print(f"{i:02d}: {geom.coords} norm(f)={norm_forces:.6f} "
        print(f"\tnorm(f)={norm_forces:.6f} "
              f"norm(ff)={norm_free_forces:.6f} dR={diff:.6f}"
        )
        if norm_free_forces < alpha:
            # Increase radius
            R += dR
            continue

        free_coords = geom.coords[free_inds]
        ref_force = forces[ref_ind]
        _ = free_forces - ref_force*(free_coords-free_ref_coords)/(constr_coord-ref_coord)
        # dir_ = _ / np.linalg.norm(_)
        # step = 0.01 * dir_
        step = _
        new_coords = geom.coords.copy()
        new_coords[free_inds] += step
        geom.coords = new_coords
    print(ref_coords)

    pot = geom.calculator
    pot.plot()
    c = np.array(c)
    rs = np.array(rs)
    pot.ax.plot(c[:,0], c[:,1], "o-")
    for rad in rs:
        circle = Circle(ref_coords[:2], radius=rad, fill=False)
        pot.ax.add_artist(circle)
    for i, xy in enumerate(c[:,:2]):
        pot.ax.annotate(i, xy)
    plt.show()


if __name__ == "__main__":
    run()
