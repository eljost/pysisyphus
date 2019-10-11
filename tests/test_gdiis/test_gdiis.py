#!/usr/bin/env python3

# [1] https://doi.org/10.1016/S0022-2860(84)87198-7
#     Pulay, 1984
# [2] https://pubs.rsc.org/en/content/articlehtml/2002/cp/b108658h
#     Farkas, Schlegel, 2002


import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot


def diis(err_vecs, max_vecs=5):
    # https://github.com/psi4/psi4numpy/blob/d8bd75a3f004728953931fb485fbc53ef8e16078/Coupled-Cluster/Spin_Orbitals/CCSD/CCSD_DIIS.py
    use_vecs = np.array(err_vecs[::-1][:max_vecs])

    # Scale error vectors so the smallest norm is 1
    norms = np.linalg.norm(use_vecs, axis=1)
    use_vecs /= norms.min()
    norms_ = np.linalg.norm(use_vecs, axis=1)

    used = len(use_vecs)

    A = np.einsum("ij,kj->ik", use_vecs, use_vecs)
    coeffs = np.linalg.solve(A, np.ones(used))
    # Scale coeffs so that their sum equals 1
    coeffs /= np.sum(coeffs)

    return coeffs, used


def test_sd_gdiis():
    geom = AnaPot.get_geom((0.4333, 3.14286, 0.))

    trust = 0.25

    cs = list()
    dcs = list()
    fs = list()
    for i in range(50):
        forces = geom.forces

        cs.append(geom.coords)
        fs.append(forces)

        forces_norm = np.linalg.norm(forces)
        print(f"{i:02d}: {forces_norm:.6f}")
        if forces_norm < 1e-3:
            print("Converged!")
            break

        if len(fs) > 1:
            coeffs, used = diis(fs, max_vecs=2)
            # Inter-/extrapolate coordinates and forces
            coords_ = np.sum(np.array(coeffs)[:,None] * cs[::-1][:used], axis=0)
            forces = np.sum(np.array(coeffs)[:,None] * fs[::-1][:used], axis=0)
            dcs.append(coords_)
            geom.coords = coords_
            forces_norm = np.linalg.norm(forces)

        step = forces
        if forces_norm > trust:
            step = trust * forces / forces_norm
        new_coords = geom.coords + step
        geom.coords = new_coords

    cs = np.array(cs)
    dcs = np.array(dcs)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*cs.T[:2], "o-")
    ax.plot(*dcs.T[:2], "o-")
    plt.show()


if __name__ == "__main__":
    test_sd_gdiis()
