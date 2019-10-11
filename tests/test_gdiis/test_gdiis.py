#!/usr/bin/env python3

# [1] https://doi.org/10.1016/S0022-2860(84)87198-7
#     Pulay, 1984
# [2] https://pubs.rsc.org/en/content/articlehtml/2002/cp/b108658h
#     Farkas, Schlegel, 2002


import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.gdiis import gdiis


def get_step(H, forces, trust):
    step = np.linalg.inv(H) @ forces
    norm = np.linalg.norm(step)
    if norm > trust:
        step = trust * step / norm
    return step


def test_sd_gdiis():
    geom = AnaPot.get_geom((0.4333, 3.14286, 0.))

    trust = 0.3
    H = np.eye(geom.coords.size)

    cs = list()
    dcs = list()
    fs = list()
    pairs = list()
    for i in range(50):
        forces = geom.forces

        cs.append(geom.coords)
        fs.append(forces)

        forces_norm = np.linalg.norm(forces)
        print(f"{i:02d}: {forces_norm:.6f}")
        if forces_norm < 1e-3:
            print("Converged!")
            break

        # Calculate reference step
        step = get_step(H, forces, trust)

        if len(fs) > 1:
            gdiis_kwargs = {
                "coords": cs,
                "forces": fs,
                "ref_step": step,
            }
            res = gdiis(fs, **gdiis_kwargs)
            if res:
                # Inter-/extrapolate coordinates and forces
                coords_ = res.coords
                forces = res.forces
                dcs.append(coords_)
                pairs.append((geom.coords.copy(), coords_, geom.coords + step))
                geom.coords = coords_
                forces_norm = np.linalg.norm(forces)
                # Get new step from DIIS coordinates & forces
                step = get_step(H, forces, trust)

        new_coords = geom.coords + step
        geom.coords = new_coords
    # return

    cs = np.array(cs)
    dcs = np.array(dcs)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*cs.T[:2], "o-")
    for i, cs_ in enumerate(cs):
        ax.annotate(f"{i:02d}", cs_[:2])

    ax.plot(*dcs.T[:2], "o-")
    for oc, dc, rs in pairs:
        line = mlines.Line2D((oc[0], dc[0]), (oc[1], dc[1]), ls="--", color="k")
        ax.add_artist(line)
        line = mlines.Line2D((oc[0], rs[0]), (oc[1], rs[1]), ls="--", color="r")
        ax.add_artist(line)
    plt.show()


if __name__ == "__main__":
    test_sd_gdiis()
