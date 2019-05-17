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
    from pysisyphus.calculators.AnaPot import AnaPot

    ref_coords = np.array((0.05, 0, 0))
    geom = CerjanMiller.get_geom(ref_coords)
    # geom = AnaPot.get_geom((-1, 1, 0))

    ref_ind = 1
    # ref_coords = geom.coords.copy()
    free_inds = np.arange(geom.coords.size) != ref_ind
    free_ref_coords = ref_coords[free_inds]
    ref_coord = ref_coords[ref_ind]

    def get_constr_coord(coords, radius):
        free_coords = coords[free_inds]
        _ = max(0, radius**2 - radius*np.sum((free_coords - free_ref_coords)**2))
        sqrt = np.sqrt(_)
        return ref_coord + sqrt

    c = list()
    dR = 0.01
    alpha = 0.01
    R = dR
    rs = [R, ]
    # import pdb; pdb.set_trace()
    for i in range(150):
        c.append(geom.coords.copy())
        print(f"Cycle {i:02d}, R={R:.04f}")
        constr_coord = get_constr_coord(geom.coords, R)
        # Set constrained coordinate
        new_coords = geom.coords.copy()
        new_coords[ref_ind] = constr_coord
        geom.coords = new_coords
        diff = R - np.linalg.norm(geom.coords - ref_coords)
        print("\t coords", geom.coords)
        geom.clear()

        # Optimize remaining coordinates
        forces = geom.forces
        print("\t forces", forces)
        norm_forces = np.linalg.norm(forces)
        # free_forces = forces[free_inds]
        # norm_free_forces = np.linalg.norm(free_forces)
        # print(f"{i:02d}: {geom.coords} norm(f)={norm_forces:.6f} "

        # free_coords = geom.coords[free_inds]
        ref_force = forces[ref_ind]
        # _ = free_forces - ref_force*(free_coords-free_ref_coords)/(constr_coord-ref_coord)
        forces_mod = forces - ref_force*(geom.coords-ref_coords)/(constr_coord-ref_coord)
        norm_forces_mod = np.linalg.norm(forces_mod)

        print(f"\tnorm(f)={norm_forces:.6f} "
              f"norm(fm)={norm_forces_mod:.6f} dR={diff:.6f}"
        )
        if norm_forces_mod < alpha:
            # Increase radius
            R += dR
            rs.append(R)
            continue

        dir_ = forces_mod / np.linalg.norm(forces_mod)
        step = 0.01 * dir_
        # step = _
        new_coords = geom.coords.copy() + step
        # new_coords[free_inds] += step
        geom.coords = new_coords
    print(ref_coords)

    pot = geom.calculator
    pot.plot()
    c = np.array(c)
    rs = np.array(rs)
    pot.ax.plot(c[:,0], c[:,1], "o-")
    pot.ax.scatter(ref_coords[0], ref_coords[1], s=20, c="r")
    for rad in rs:
        circle = Circle(ref_coords[:2], radius=rad, fill=False)
        pot.ax.add_artist(circle)
    for i, xy in enumerate(c[:,:2]):
        pot.ax.annotate(i, xy)
    # pot.ax.set_xlim(-0.2, 0.2)
    # pot.ax.set_ylim(-0.2, 0.2)
    plt.show()

    # from pysisyphus.calculators.AnaPot import AnaPot
    # g = AnaPot.get_geom((-1, 1, 0))
    # g.calculator.plot()
    # plt.show()


if __name__ == "__main__":
    run()
