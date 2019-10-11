#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot


def diis(err_vecs, max_vecs=5):
    # https://github.com/psi4/psi4numpy/blob/d8bd75a3f004728953931fb485fbc53ef8e16078/Coupled-Cluster/Spin_Orbitals/CCSD/CCSD_DIIS.py
    use_vecs = np.array(err_vecs[::-1][:max_vecs])
    used = len(use_vecs)
    print(f"using last {len(use_vecs)} error vectors")

    size_ = len(use_vecs)+1
    A = np.ones((size_, size_))
    A[-1,-1] = 0
    _ = np.einsum("ij,jk->ik", use_vecs, use_vecs.T)
    A[:-1,:-1] = _
    sol = np.zeros(size_)
    sol[-1] = 1
    *coeffs, lagrange = np.linalg.solve(A, sol)
    print(A, coeffs, lagrange)
    return coeffs, used


def test_sd_gdiis():
    geom = AnaPot.get_geom((0.4333, 3.14286, 0.))

    trust = 0.25

    cs = list()
    dcs = list()
    gs = list()
    for i in range(50):
        forces = geom.forces

        cs.append(geom.coords)
        gs.append(-forces)

        forces_norm = np.linalg.norm(forces)
        print(f"{i:02d}: {forces_norm:.6f}")
        if forces_norm < 1e-3:
            print("Converged!")
            break

        if len(gs) > 1:
            coeffs, used = diis(gs, max_vecs=2)
            # Inter-/extrapolate coordinates and gradient
            coords_ = np.sum(np.array(coeffs)[:,None] * cs[::-1][:used], axis=0)
            forces = -np.sum(np.array(coeffs)[:,None] * gs[::-1][:used], axis=0)
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

    pass


if __name__ == "__main__":
    test_sd_gdiis()
